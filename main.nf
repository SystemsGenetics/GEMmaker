#!/usr/bin/env nextflow
/*
================================================================================
                         GENmaker
================================================================================
 GEMmaker Analysis Pipeline.
 #### Homepage / Documentation
 https://github.com/SystemsGenetics/gemmaker
 https://gemmaker.readthedocs.io/en/latest/
--------------------------------------------------------------------------------
*/

import java.nio.channels.FileLock
import java.nio.channels.FileChannel
import java.nio.channels.OverlappingFileLockException

log.info Headers.nf_core(workflow, params.monochrome_logs)

////////////////////////////////////////////////////
/* --               PRINT HELP                 -- */
////////////////////////////////////////////////////+
def json_schema = "$projectDir/nextflow_schema.json"
if (params.help) {
    def command = "nextflow run systemsgenetics/gemmaker -profile singularity --pipeline kallisto --kallisto_index_path Arabidopsis_thaliana.TAIR10.kallisto.indexed --sras SRAs.txt"
    log.info NfcoreSchema.params_help(workflow, params, json_schema, command)
    exit 0
}

////////////////////////////////////////////////////
/* --         VALIDATE PARAMETERS              -- */
////////////////////////////////////////////////////+
if (params.validate_params) {
    NfcoreSchema.validateParameters(params, json_schema, log)
}

////////////////////////////////////////////////////
/* --     Collect configuration parameters     -- */
////////////////////////////////////////////////////

println """\
Workflow Information:
---------------------
  Project Directory:          ${workflow.projectDir}
  Work Directory:             ${workflow.workDir}
  Launch Directory:           ${workflow.launchDir}
  Config Files:               ${workflow.configFiles}
  Container Engine:           ${workflow.containerEngine}
  Profile(s):                 ${workflow.profile}


Samples:
--------
  Remote fastq list path:     ${params.sras}
  Local sample glob:          ${params.input}
  Skip samples file:          ${params.skip_samples}


Reports
-------
  Report directory:           ${params.outdir}/reports


Quantification:
---------------
  Tool:                       ${params.pipeline}"""

// Indicates which tool the user selected.
hisat2_enable = false
kallisto_enable = false
salmon_enable = false
selected_tool = 0

// Print out details per the selected tool.
if (params.pipeline.equals('hisat2')) {
  hisat2_enable = true
  selected_tool = 0
  println """\
  Hisat2 Index Base Name:     ${params.hisat2_base_name}
  Hisat2 GTF File:            ${params.hisat2_gtf_file}
  Hisat2 Index Directory:     ${params.hisat2_index_dir}

  Trimmomatic Parameters:
     clip path:               ${params.trimmomatic_clip_file}
     MINLEN:                  ${params.trimmomatic_MINLEN}
     SLIDINGWINDOW:           ${params.trimmomatic_SLIDINGWINDOW}
     LEADING:                 ${params.trimmomatic_LEADING}
     TRAILING:                ${params.trimmomatic_TRAILING}
  """
}
if (params.pipeline.equals('kallisto')) {
  kallisto_enable = true
  selected_tool = 1
  println """\
  Kallisto Index File:        ${params.kallisto_index_path}
  """
}
if (params.pipeline.equals('salmon')) {
  salmon_enable = true
  selected_tool = 2
  println """\
  Salmon Index Directory:     ${params.salmon_index_path}
  """
}

if (!(hisat2_enable || kallisto_enable || salmon_enable)) {
  error "Error: You must select a valid quantification tool using the '--pipeline' parameter. Currently valid options are 'salmon', 'kallisto', or 'hisat2'"
}

// Create the directories we'll use for running
// batches
file("${workflow.workDir}/GEMmaker").mkdir()
file("${workflow.workDir}/GEMmaker/stage").mkdir()
file("${workflow.workDir}/GEMmaker/process").mkdir()
file("${workflow.workDir}/GEMmaker/done").mkdir()

// Make sure that once GEMmaker runs that the user doesn't try to change
// quantification tools half way through.
method_lock_file = file("${workflow.workDir}/GEMmaker/method")
if (!method_lock_file.exists()) {
    method_lock_file << "${params.pipeline}"
}
else {
  reader = method_lock_file.newReader()
  active_method = reader.readLine()
  reader.close()
  if (!active_method.equals(params.pipeline)) {
    error "Error: previously, GEMmaker was set to run using the '${active_method}' tool, but it looks as though the configuration has changed to use the '${params.pipeline}' tool. GEMmaker only supports use of one tool at a time. If you would like to change the quantification tool please re-run GEMmaker in a new directory or remove the `work` and `results` directories prior to restarting GEMmaker to clear out unwanted results."
  }
}

/**
 * Check to make sure that required reference files exist
 */
// If Hisat2 was selected:
if (hisat2_enable) {
  gtfFile = file("${params.hisat2_gtf_file}")

  if (gtfFile.isEmpty()) {
    error "Error: GTF reference file for Hisat2 does not exist or is empty! Please Check that you have the proper references, that they are placed in the reference directory, and they are named properly.\
    \nGEMmaker is missing the following file: '${params.hisat2_gtf_file}' "
  }

  hisat2_index_dir = file("${params.hisat2_index_dir}")

  if (!hisat2_index_dir.isDirectory()) {
    error "Error: hisat2 Index Directory does not exist or is empty! Please Check that you have the proper references, that they are placed in the reference directory, and they are named properly.\
    \nGEMmaker is missing the following file: '${params.hisat2_index_dir}'"
  }
}

// If Kallisto was selected
if (kallisto_enable) {
  kallisto_index_file = file("${params.kallisto_index_path}")

  if (kallisto_index_file.isEmpty()) {
    error "Error: Kallisto Index File does not exist or is empty! Please Check that you have the proper references, that they are placed in the reference directory, and they are named properly.\
    \nGEMmaker is missing the following file: '${params.kallisto_index_path}'"
  }
}

// If Salmon was selected
if (salmon_enable) {
  salmon_index_dir = file("${params.salmon_index_path}")

  if (!salmon_index_dir.isDirectory()) {
    error "Error: Salmon Index Directory does not exist or is empty! Please Check that you have the proper references, that they are placed in the reference directory, and they are named properly.\
    \nGEMmaker is missing the following file: '${params.salmon_index_path}'"
  }
}

/**
 * Check that other input files/directories exist
 */
 if (hisat2_enable) {
     clip_path = file("${params.trimmomatic_clip_file}")
     if (!clip_path.exists()) {
       error "Error: The Trimmomatic clip file cannot be found at '${params.trimmomatic_clip_file}'."
     }
 }

if (params.sras) {
    sample_file = file("${params.sras}")
    if (!sample_file.exists()) {
        error "Error: The NCBI download sample file does not exists at '${params.sras}'. This file must be provided. If you are not downloading samples from NCBI SRA the file must exist but can be left empty."
    }
}

if (params.skip_samples) {
    skip_file = file("${params.skip_samples}")
    if (!skip_file.exists()) {
       error "Error: The file containing samples to skip does not exists at '${params.skip_samples}'."
   }
}

failed_report_template = file("${params.failed_run_report_template}")
if (!failed_report_template.exists()) {
   error "Error: The failed run report template cannot be found at '${params.failed_run_report_template}'. This file comes with GEMmaker and is required."
}

/**
 * Create value channels that can be reused
 */
if (hisat2_enable) {
    HISAT2_INDEXES = Channel.fromPath("${params.hisat2_index_dir}/*").collect()
    FASTA_ADAPTER = Channel.fromPath("${params.trimmomatic_clip_file}").collect()
    GTF_FILE = Channel.fromPath("${params.hisat2_gtf_file}").collect()
}
else {
    Channel.empty().set { HISAT2_INDEXES }
    Channel.empty().set { FASTA_ADAPTER }
    Channel.empty().set { GTF_FILE }
}
if (kallisto_enable) {
    KALLISTO_INDEX = Channel.fromPath("${params.kallisto_index_path}").collect()
}
else {
    Channel.empty().set { KALLISTO_INDEX }
}
if (salmon_enable) {
    SALMON_INDEXES = Channel.fromPath("${params.salmon_index_path}").collect()
}
else {
    Channel.empty().set { SALMON_INDEXES }
}
FAILED_RUN_TEMPLATE = Channel.fromPath("${params.failed_run_report_template}").collect()
MULTIQC_CONFIG = Channel.fromPath("${params.multiqc_config_file}").collect()
MULTIQC_CUSTOM_LOGO = Channel.fromPath("${params.multiqc_custom_logo}").collect()

if (params.skip_samples) {
  SKIP_SAMPLES_FILE = Channel.fromPath("${params.skip_samples}")
}
else {
    Channel.empty().set { SKIP_SAMPLES_FILE }
}

/**
 * Local Sample Input.
 * This checks the folder that the user has given
 */
if (!params.input) {
  Channel.empty().set { LOCAL_SAMPLE_FILES_FOR_STAGING }
  Channel.empty().set { LOCAL_SAMPLE_FILES_FOR_JOIN }
}
else {
  Channel.fromFilePairs( "${params.input}", size: -1 )
    .set { LOCAL_SAMPLE_FILES_FOR_STAGING }
  Channel.fromFilePairs( "${params.input}", size: -1 )
    .set { LOCAL_SAMPLE_FILES_FOR_JOIN }
}

/**
 * Remote fastq_run_id Input.
 */
if (params.sras == "") {
  Channel.empty().set { SRR_FILE }
}
else {
  Channel.fromPath("${params.sras}").set { SRR_FILE }
}



/**
 * Make sure that at least one output format is enabled.
 */
if (hisat2_enable == true && params.hisat2_keep_counts == false &&
    params.hisat2_keep_fpkm == false && params.hisat2_keep_tpm == false) {
  error "Error: at least one output format (raw, fpkm, tpm) must be enabled for hisat2"
}

if (hisat2_enable == false && params.hisat2_keep_counts == false && params.hisat2_keep_tpm == false) {
  error "Error: at least one output format (raw, tpm) must be enabled for kallisto / salmon"
}

println """\

Published Results:
---------------
  Output Dir:                 ${params.outdir}/GEMs """

// For the create_gem process we need to know if the FKPM, TMP and raw counts
// should be published.
publish_fpkm = false
publish_tpm = false
publish_raw = false
publish_gem = false
if (hisat2_enable && params.hisat2_keep_fpkm) {
    publish_fpkm = true
    println """  FPKM counts:                Yes"""
}
if ((hisat2_enable && params.hisat2_keep_counts) ||
    (salmon_enable && params.salmon_keep_counts) ||
    (kallisto_enable && params.kallisto_keep_counts)) {
    publish_raw = true
    println """  Raw counts:                 Yes"""
}
if ((hisat2_enable && params.hisat2_keep_tpm) ||
    (salmon_enable && params.salmon_keep_tpm) ||
    (kallisto_enable && params.kallisto_keep_tpm)) {
    publish_tpm = true
    println """  TPM counts:                 Yes"""
}

if ((hisat2_enable && params.hisat2_keep_gem) ||
    (salmon_enable && params.salmon_keep_gem) ||
    (kallisto_enable && params.kallisto_keep_gem)) {
    publish_gem = true
    println """  GEM file:                   Yes"""
}

// Add a few lines after the header.
println """\

"""


/**
 * Set the pattern for publishing downloaded FASTQ files
 */
publish_pattern_fastq_dump = params.keep_retrieved_fastq
  ? "{*.fastq}"
  : "{none}"

/**
 * Set the pattern for publishing trimmed FASTQ files
 */
publish_pattern_trimmomatic = params.trimmomatic_keep_trimmed_fastq
  ? "{*.trim.log,*_trim.fastq}"
  : "{*.trim.log}"

/**
 * Set the pattern for publishing BAM files
 */
publish_pattern_samtools_sort = params.hisat2_keep_bam
  ? "{*.log,*.bam}"
  : "{*.log}"

publish_pattern_samtools_index = params.hisat2_keep_bam
  ? "{*.log,*.bam.bai}"
  : "{*.log}"

/**
 * Set the pattern for publishing Kallisto GA files
 */
publish_pattern_Kallisto_GA = params.kallisto_keep_data
  ? "{*.ga,*.log}"
  : "{*.log}"

/**
 * Set the pattern for publishing Salmon GA files
 * Publishes only log file used by multiqc if false
 */
publish_pattern_Salmon_GA = params.salmon_keep_data
  ? "{*.ga}"
  : "{*.ga/aux_info/meta_info.json,*.ga/libParams/flenDist.txt}"

/**
 * Set the pattern for publishing STRINGTIE GA and GTF files
 */
publish_pattern_stringtie_gtf_and_ga = params.hisat2_keep_data
  ? "{*.ga, *.gtf}"
  : "{none}"



/**
 * Retrieves metadata for all of the remote samples
 * and maps SRA runs to SRA experiments.
 */
process retrieve_sra_metadata {
  publishDir params.outdir, mode: params.publish_dir_mode, pattern: "failed_runs.metadata.txt"
  label "retrieve_sra_metadata"

  input:
    file srr_file from SRR_FILE
    file skip_samples from SKIP_SAMPLES_FILE

  output:
    stdout REMOTE_SAMPLES_LIST
    file "failed_runs.metadata.txt" into METADATA_FAILED_RUNS

  script:
  if (skip_samples) {
      skip_arg = "--skip_file ${skip_samples}"
  }
  """
  >&2 echo "#TRACE n_remote_run_ids=`cat ${srr_file} | wc -l`"

  retrieve_sra_metadata.py \
      --run_id_file ${srr_file} \
      --meta_dir ${workflow.workDir}/GEMmaker \
      ${skip_arg}
  """
}



/**
 * Splits the SRR2XRX mapping file
 */

// First create a list of the remote and local samples
REMOTE_SAMPLES_LIST
  .splitCsv()
  .groupTuple(by: 1)
  .map { [it[1], it[0].toString().replaceAll(/[\[\]\'\,]/,''), 'remote'] }
  .set{REMOTE_SAMPLES_FOR_STAGING}

LOCAL_SAMPLE_FILES_FOR_STAGING
  .map{ [it[0], it[1], 'local'] }
  .set{LOCAL_SAMPLES_FOR_STAGING}

ALL_SAMPLES = REMOTE_SAMPLES_FOR_STAGING
  .mix(LOCAL_SAMPLES_FOR_STAGING)


// Channels to bootstrap post-processing of
// sample results if a resume is performed when
// all samples have completed.
MULTIQC_BOOTSTRAP = Channel.create()
CREATE_GEM_BOOTSTRAP = Channel.create()

// Remove any lock file that might be leftover from a previous run
lockfile = file("${workflow.workDir}/GEMmaker/gemmaker.lock")
if (lockfile.exists()) {
    lockfile.delete()
}

// Clean up any files left over from a previous run by moving them
// back to the stage directory.
existing_files = file('work/GEMmaker/process/*')
for (existing_file in existing_files) {
  existing_file.moveTo('work/GEMmaker/stage')
}

// Check to see if we have any files left in the
// stage directory. If so we need to keep processing
// samples
staged_files = file('work/GEMmaker/stage/*')

// If a user added a sample to skip after a failed run, then
// we want to remove it from the stage folder.
if (params.skip_samples) {
    skip_file = file("${params.skip_samples}")
    skip_file.eachLine { line ->
        skip_sample = file('work/GEMmaker/stage/' + line.trim() + '.sample.csv')
        if (skip_sample.exists()) {
            skip_sample.delete()
        }
    }
}


if (staged_files.size() == 0) {
  // If there are no staged files then the workflow will
  // end because it only proceeds when there are samples
  // in the processed directory.  However suppose the workflow
  // fails on multiqc and needs to be resumed.  The
  // following bootstraps the post-processsing portion of
  // the workflow
  MULTIQC_BOOTSTRAP.bind(1)
  CREATE_GEM_BOOTSTRAP.bind(1)
}


/**
 * Writes the batch files and stores them in the
 * stage directory.
 */
process write_stage_files {
  tag { sample_id }
  label "local"

  input:
    set val(sample_id), val(run_files_or_ids), val(sample_type) from ALL_SAMPLES

  output:
    val (1) into SAMPLES_READY_SIGNAL

  exec:
    // Get any samples to skip
    skip_samples = []
    if (params.skip_samples) {
        skip_file = file("${params.skip_samples}")
        skip_file.eachLine { line ->
            skip_samples << line.trim()
        }
    }

    // Only stage files that should not be skipped.
    if (skip_samples.intersect([sample_id]) == []) {
      // Create a file for each samples.
      sample_file = file("${workflow.workDir}/GEMmaker/stage/" + sample_id + '.sample.csv')
      sample_file.withWriter {
        // If this is a local file.
        if (sample_type.equals('local')) {
          it.writeLine '"' + sample_id + '","' + run_files_or_ids.join('::') + '","' + sample_type + '"'
        }
        // If this is a remote file.
        else {
          it.writeLine '"' + sample_id + '","' + run_files_or_ids + '","' + sample_type + '"'
        }
      }
    }
}



// When all batch files are created we need to then
// move the first file into the process directory.
SAMPLES_READY_SIGNAL.collect().set { FIRST_SAMPLE_START_SIGNAL }

/*
 * Parse software version numbers
 */
process get_software_versions {
    publishDir "${params.outdir}/pipeline_info", mode: params.publish_dir_mode,
        saveAs: { filename ->
                      if (filename.indexOf('.csv') > 0) filename
                      else null
        }

    output:
    file 'software_versions_mqc.yaml' into ch_software_versions_yaml
    file 'software_versions.csv'

    script:
    """
    echo $workflow.manifest.version > v_pipeline.txt
    echo $workflow.nextflow.version > v_nextflow.txt
    fastqc --version > v_fastqc.txt
    multiqc --version > v_multiqc.txt
    kallisto version > v_kallisto.txt
    hisat2 --version | head -n 1 > v_hisat2.txt
    salmon --version > v_salmon.txt
    python --version > v_python.txt
    samtools version | head -n 1 > v_samtools.txt
    fastq-dump --version > v_fastq_dump.txt
    stringtie --version > v_stringtie.txt
    trimmomatic -version > v_trimmomatic.txt
    pip freeze | grep pandas > v_pandas.txt
    scrape_software_versions.py &> software_versions_mqc.yaml
    """
}


/**
 * Moves the first set of sample files into the process directory.
 */
process start_first_batch {
  label "local"
  cache false

  input:
    val signal from FIRST_SAMPLE_START_SIGNAL

  exec:
    // Move the first set of sample file into the processing directory
    // so that we jumpstart the workflow.
    sample_files = file("${workflow.workDir}/GEMmaker/stage/*.sample.csv")
    start_samples = sample_files.sort().take(params.max_cpus)
    if (sample_files.size() > 0 ) {
      for (sample in start_samples) {
        sample.moveTo("${workflow.workDir}/GEMmaker/process")
      }
   }
   // If there are no staged files then we need to
   // close out the channels so we don't hang.
   else {
      NEXT_SAMPLE.close()
      NEXT_SAMPLE_SIGNAL.close()
      HISAT2_SAMPLE_COMPLETE_SIGNAL.close()
      KALLISTO_SAMPLE_COMPLETE_SIGNAL.close()
      SALMON_SAMPLE_COMPLETE_SIGNAL.close()
      SAMPLE_COMPLETE_SIGNAL.close()
      MULTIQC_BOOTSTRAP.close()
      CREATE_GEM_BOOTSTRAP.close()
      SKIP_DUMP_SAMPLE.close()
      SKIP_DOWNLOAD_SAMPLE.close()
      println "There are no staged samples.  Moving on to post-processing"
   }
}



// Create the channel that will watch the process directory
// for new files. When a new sample file is added
// it will be read it and sent it through the workflow.
NEXT_SAMPLE = Channel
   .watchPath("${workflow.workDir}/GEMmaker/process")



/**
 * Opens the sample file and prints it's contents to
 * STDOUT so that the samples can be caught in a new
 * channel and start processing.
 */
process read_sample_file {
  tag { sample_file }
  label "local"
  cache false

  input:
    file(sample_file) from NEXT_SAMPLE

  output:
    stdout SAMPLE_FILE_CONTENTS

  script:
  """
  cat ${sample_file}
  """
}



// Split our sample file contents into two different
// channels, one for remote samples and another for local.
LOCAL_SAMPLES = Channel.create()
REMOTE_SAMPLES = Channel.create()
SAMPLE_FILE_CONTENTS
  .splitCsv(quote: '"')
  .choice(LOCAL_SAMPLES, REMOTE_SAMPLES) { a -> a[2] =~ /local/ ? 0 : 1 }

// Split our list of local samples into two pathways, onefor
// FastQC analysis and the other for read counting.  We don't
// do this for remote samples because they need downloading
// first.
LOCAL_SAMPLES
  .map {[it[0], 'hi']}
  .mix(LOCAL_SAMPLE_FILES_FOR_JOIN)
  .groupTuple(size: 2)
  .map {[it[0], it[1][0]]}
  .into {LOCAL_SAMPLES_FOR_FASTQC_1; LOCAL_SAMPLES_FOR_COUNTING}



// Create the channels needed for signalling when
// samples are completed.
HISAT2_SAMPLE_COMPLETE_SIGNAL = Channel.create()
KALLISTO_SAMPLE_COMPLETE_SIGNAL = Channel.create()
SALMON_SAMPLE_COMPLETE_SIGNAL = Channel.create()
SKIP_DOWNLOAD_SAMPLE = Channel.create()
SKIP_DUMP_SAMPLE = Channel.create()

// Create the channel that will collate all the signals
// and release a signal when the sample is complete
SAMPLE_COMPLETE_SIGNAL = Channel.create()
SAMPLE_COMPLETE_SIGNAL
  .mix(HISAT2_SAMPLE_COMPLETE_SIGNAL, KALLISTO_SAMPLE_COMPLETE_SIGNAL, SALMON_SAMPLE_COMPLETE_SIGNAL,
      SKIP_DUMP_SAMPLE.map {it[0]}, SKIP_DOWNLOAD_SAMPLE.map {it[0]})
  .into { NEXT_SAMPLE_SIGNAL; MULTIQC_READY_SIGNAL; CREATE_GEM_READY_SIGNAL }



/**
 * Handles the end of a sample by moving a new sample
 * file into the process directory which triggers
 * the NEXT_SAMPLE.watchPath channel.
 */
process next_sample {
  tag { sample_id }
  label "local"
  cache false

  input:
    val sample_id from NEXT_SAMPLE_SIGNAL

  exec:
    // Use a file lock to prevent a race condition for grabbing the next sample.
    FileChannel channel = null
    FileLock lock = null
    success = false

    try {
      attempts = 0
      channel = new RandomAccessFile("${workflow.workDir}/GEMmaker/gemmaker.lock", "rw").getChannel()
      // We will try for a release lock for about 1 hour. It's not clear
      // what the optimal wait time should be.
      while (!lock)  {
        // 60 attempts with a sleep of 1 minute == 1 hour.
        if (attempts < 60) {
          try {
            lock = channel.lock()
          }
          catch (OverlappingFileLockException e) {
            // Do nothing, let's try a few more times....
          }
          if (!lock) {
            println "Waiting on lock. After sample, " + sample_id + ", attempt " + attempts + "..."
            // Sleep for 1 minute
            sleep 60000
            attempts = attempts + 1
          }
        }
        else {
          throw new Exception("Cannot obtain lock to proceed to next sample after 3 attempts")
        }
      }
      sample_file = file("${workflow.workDir}/GEMmaker/process/" + sample_id + '.sample.csv')
      sample_file.moveTo("${workflow.workDir}/GEMmaker/done")

      // Move the next sample file into the processing directory
      // which will trigger the start of the next sample.
      staged_files = file("${workflow.workDir}/GEMmaker/stage/*")
      if (staged_files.size() > 0) {
        staged_files.first().moveTo("${workflow.workDir}/GEMmaker/process")
      }
      else {
        processing_files = file("${workflow.workDir}/GEMmaker/process/*.sample.csv")
        if (processing_files.size() == 0) {
          NEXT_SAMPLE.close()
          NEXT_SAMPLE_SIGNAL.close()
          HISAT2_SAMPLE_COMPLETE_SIGNAL.close()
          KALLISTO_SAMPLE_COMPLETE_SIGNAL.close()
          SALMON_SAMPLE_COMPLETE_SIGNAL.close()
          SAMPLE_COMPLETE_SIGNAL.close()
          MULTIQC_BOOTSTRAP.close()
          CREATE_GEM_BOOTSTRAP.close()
          SKIP_DUMP_SAMPLE.close()
          SKIP_DOWNLOAD_SAMPLE.close()
        }
      }
      success = true
    }
    catch (Exception e) {
      println "Error: " + e.getMessage()
    }
    finally {
      // Release the lock file and close the file if they were opened.
      if (lock && lock.isValid()) {
        lock.release()
      }
      // Close the channel.
      if (channel) {
        channel.close()
      }
      lock = null;
      // Re-throw exception to terminate the workflow if there was no success.
      if (!success) {
        throw new Exception("Could not move to the next sample.")
      }
    }
}



/**
 * Downloads SRA files from NCBI using the SRA Toolkit.
 */
process download_runs {
  publishDir params.outdir, mode: params.publish_dir_mode, pattern: '*.failed_runs.download.txt', saveAs: { "Samples/${sample_id}/${it}" }

  tag { sample_id }
  label "download_runs"

  input:
    set val(sample_id), val(run_ids), val(type) from REMOTE_SAMPLES

  output:
    set val(sample_id), file("*.sra") optional true into SRA_TO_EXTRACT
    set val(sample_id), file("*.sra") optional true into SRA_TO_CLEAN
    set val(sample_id), file("sample_failed") optional true into SKIP_DOWNLOAD_SAMPLE
    set val(sample_id), file('*.failed_runs.download.txt') into DOWNLOAD_FAILED_RUNS

  script:
  """
  echo "#TRACE n_remote_run_ids=${run_ids.tokenize(' ').size()}"

  retrieve_sra.py --sample ${sample_id} --run_ids ${run_ids} --akey \$ASPERA_KEY
  """
}



/**
 * Extracts FASTQ files from downloaded SRA files.
 */
process fastq_dump {
  publishDir params.outdir, mode: params.publish_dir_mode, pattern: publish_pattern_fastq_dump, saveAs: { "Samples/${sample_id}/${it}" }
  publishDir params.outdir, mode: params.publish_dir_mode, pattern: '*.failed_runs.fastq-dump.txt', saveAs: { "Samples/${sample_id}/${it}" }
  tag { sample_id }
  label "fastq_dump"

  input:
    set val(sample_id), file(sra_files) from SRA_TO_EXTRACT

  output:
    set val(sample_id), file("*.fastq") optional true into DOWNLOADED_FASTQ_FOR_MERGING
    set val(sample_id), file("*.fastq") optional true into DOWNLOADED_FASTQ_FOR_CLEANING
    set val(sample_id), file("sample_failed") optional true into SKIP_DUMP_SAMPLE
    set val(sample_id), val(1) into CLEAN_SRA_SIGNAL
    set val(sample_id), file('*.failed_runs.fastq-dump.txt') into FASTQ_DUMP_FAILED_RUNS

  script:
  """
  echo "#TRACE sample_id=${sample_id}"
  echo "#TRACE sra_bytes=`stat -Lc '%s' *.sra | awk '{sum += \$1} END {print sum}'`"

  sra2fastq.py --sample ${sample_id} --sra_files ${sra_files}
  """
}



/**
 * This process merges the fastq files based on their sample_id number.
 */
process fastq_merge {
  tag { sample_id }

  input:
    set val(sample_id), file(fastq_files) from DOWNLOADED_FASTQ_FOR_MERGING

  output:
    set val(sample_id), file("${sample_id}_?.fastq") into MERGED_SAMPLES_FOR_COUNTING
    set val(sample_id), file("${sample_id}_?.fastq") into MERGED_SAMPLES_FOR_FASTQC_1
    set val(sample_id), file("${sample_id}_?.fastq") into MERGED_FASTQ_FOR_CLEANING
    set val(sample_id), val(1) into CLEAN_FASTQ_SIGNAL

  script:
  """
  echo "#TRACE sample_id=${sample_id}"
  echo "#TRACE fastq_lines=`cat *.fastq | wc -l`"

  fastq_merge.sh ${sample_id}
  """
}



/**
 * This is where we combine samples from both local and remote sources.
 */
COMBINED_SAMPLES_FOR_FASTQC_1 = LOCAL_SAMPLES_FOR_FASTQC_1.mix(MERGED_SAMPLES_FOR_FASTQC_1)
COMBINED_SAMPLES_FOR_COUNTING = LOCAL_SAMPLES_FOR_COUNTING.mix(MERGED_SAMPLES_FOR_COUNTING)

/**
 * Performs fastqc on raw fastq files
 */
process fastqc_1 {
  publishDir "${params.outdir}/Samples/${sample_id}", mode: params.publish_dir_mode, pattern: "*_fastqc.*"
  tag { sample_id }
  label "fastqc"

  input:
    set val(sample_id), file(fastq_files) from COMBINED_SAMPLES_FOR_FASTQC_1

  output:
    set file("*_fastqc.html") , file("*_fastqc.zip") into FASTQC_1_OUTPUT
    set val(sample_id), val(1) into CLEAN_MERGED_FASTQ_FASTQC_SIGNAL

  script:
  """
  echo "#TRACE sample_id=${sample_id}"
  echo "#TRACE fastq_lines=`cat *.fastq | wc -l`"

  fastqc ${fastq_files}
  """
}



/**
 * THIS IS WHERE THE SPLIT HAPPENS FOR hisat2 vs Kallisto vs Salmon
 *
 * Information about "choice" split operator (to be deleted before final
 * GEMmaker release)
 */
HISAT2_CHANNEL = Channel.create()
KALLISTO_CHANNEL = Channel.create()
SALMON_CHANNEL = Channel.create()
COMBINED_SAMPLES_FOR_COUNTING.choice( HISAT2_CHANNEL, KALLISTO_CHANNEL, SALMON_CHANNEL) { selected_tool }



/**
 * Performs KALLISTO alignemnt of fastq files
 */
process kallisto {
  publishDir "${params.outdir}/Samples/${sample_id}", mode: params.publish_dir_mode, pattern: publish_pattern_Kallisto_GA
  tag { sample_id }
  label "multithreaded"
  label "kallisto"

  input:
    set val(sample_id), file(fastq_files) from KALLISTO_CHANNEL
    file kallisto_index from KALLISTO_INDEX

  output:
    set val(sample_id), file("*.ga") into KALLISTO_GA
    set val(sample_id), file("*.ga") into KALLISTO_GA_TO_CLEAN
    set val(sample_id), val(1) into CLEAN_MERGED_FASTQ_KALLISTO_SIGNAL
    file "*.kallisto.log" into KALLISTO_LOG

  script:
  """
  echo "#TRACE sample_id=${sample_id}"
  echo "#TRACE fastq_lines=`cat *.fastq | wc -l`"
  echo "#TRACE index_bytes=`stat -Lc '%s' ${kallisto_index}`"

  kallisto.sh \
    ${sample_id} \
    ${kallisto_index} \
    ${task.cpus} \
    "${fastq_files}"
  """
}



/**
 * Generates the final TPM and raw count files for Kallisto
 */
process kallisto_tpm {
  publishDir "${params.outdir}/Samples/${sample_id}", mode: params.publish_dir_mode
  tag { sample_id }

  input:
    set val(sample_id), file(ga_File) from KALLISTO_GA

  output:
    file "*.tpm" optional true into KALLISTO_TPM
    file "*.raw" optional true into KALLISTO_RAW
    set val(sample_id), val(1) into CLEAN_KALLISTO_GA_SIGNAL
    val sample_id into KALLISTO_SAMPLE_COMPLETE_SIGNAL

  script:
  """
  echo "#TRACE sample_id=${sample_id}"
  #echo "#TRACE ga_lines=`cat *.ga | wc -l`"

  kallisto_tpm.sh \
    ${sample_id} \
    ${params.kallisto_keep_tpm} \
    ${params.kallisto_keep_counts}
  """
}



/**
 * Performs SALMON alignemnt of fastq files
 */
process salmon {
  publishDir "${params.outdir}/Samples/${sample_id}", mode: params.publish_dir_mode, pattern: publish_pattern_Salmon_GA
  tag { sample_id }
  label "multithreaded"
  label "salmon"

  input:
    set val(sample_id), file(fastq_files) from SALMON_CHANNEL
    file salmon_index from SALMON_INDEXES

  output:
    set val(sample_id), file("*.ga") into SALMON_GA
    set val(sample_id), path("${sample_id}.Salmon.ga", type: 'dir') into SALMON_GA_LOG
    set val(sample_id), file("*.ga/quant.sf") into SALMON_GA_TO_CLEAN
    set val(sample_id), val(1) into CLEAN_MERGED_FASTQ_SALMON_SIGNAL

  script:
  """
  echo "#TRACE sample_id=${sample_id}"
  echo "#TRACE fastq_lines=`cat *.fastq | wc -l`"
  echo "#TRACE index_bytes=`stat -Lc '%s' ${salmon_index} | awk '{sum += \$1} END {print sum}'`"

  salmon.sh \
    ${sample_id} \
    ${task.cpus} \
    ${salmon_index} \
    "${fastq_files}" \
  """
}



/**
 * Generates the final TPM file for Salmon
 */
process salmon_tpm {
  publishDir "${params.outdir}/Samples/${sample_id}", mode: params.publish_dir_mode
  tag { sample_id }

  input:
    set val(sample_id), file(ga_file) from SALMON_GA

  output:
    file "*.Salmon.tpm" optional true into SALMON_TPM
    file "*.Salmon.raw" optional true into SALMON_RAW
    set val(sample_id), val(1) into CLEAN_SALMON_GA_SIGNAL
    val sample_id into SALMON_SAMPLE_COMPLETE_SIGNAL

  script:
  """
  echo "#TRACE sample_id=${sample_id}"

  salmon_tpm.sh \
    ${params.salmon_keep_tpm} \
    ${params.salmon_keep_counts} \
    ${sample_id}
  """
}



/**
 * Performs Trimmomatic on all fastq files.
 *
 * This process requires that the ILLUMINACLIP_PATH environment
 * variable be set in the trimmomatic module. This indicates
 * the path where the clipping files are stored.
 *
 * MINLEN is calculated using based on percentage of the mean
 * read length. The percentage is determined by the user in the
 * "nextflow.config" file
 */
process trimmomatic {
  publishDir "${params.outdir}/Samples/${sample_id}", mode: params.publish_dir_mode, pattern: publish_pattern_trimmomatic
  tag { sample_id }
  label "multithreaded"
  label "trimmomatic"

  input:
    set val(sample_id), file(fastq_files) from HISAT2_CHANNEL
    file fasta_adapter from FASTA_ADAPTER

  output:
    set val(sample_id), file("*_trim.fastq") into TRIMMED_SAMPLES_FOR_FASTQC
    set val(sample_id), file("*_trim.fastq") into TRIMMED_SAMPLES_FOR_HISAT2
    set val(sample_id), file("*_trim.fastq") into TRIMMED_FASTQ_FOR_CLEANING
    set val(sample_id), file("*.trim.log") into TRIMMED_SAMPLE_LOG

  script:
  """
  echo "#TRACE sample_id=${sample_id}"
  echo "#TRACE n_cpus=${task.cpus}"
  echo "#TRACE minlen=${params.trimmomatic_MINLEN}"
  echo "#TRACE leading=${params.trimmomatic_LEADING}"
  echo "#TRACE trailing=${params.trimmomatic_TRAILING}"
  echo "#TRACE slidingwindow=${params.trimmomatic_SLIDINGWINDOW}"
  echo "#TRACE fasta_lines=`cat ${fasta_adapter} | wc -l`"
  echo "#TRACE fastq_lines=`cat *.fastq | wc -l`"

  trimmomatic.sh \
    ${sample_id} \
    ${params.trimmomatic_MINLEN} \
    ${task.cpus} \
    ${fasta_adapter} \
    ${params.trimmomatic_LEADING} \
    ${params.trimmomatic_TRAILING} \
    ${params.trimmomatic_SLIDINGWINDOW} \
    "${fastq_files}"
  """
}



/**
 * Performs fastqc on fastq files post trimmomatic
 * Files are stored to an independent folder
 */
process fastqc_2 {
  publishDir "${params.outdir}/Samples/${sample_id}", mode: params.publish_dir_mode, pattern: "*_fastqc.*"
  tag { sample_id }
  label "fastqc"

  input:
    set val(sample_id), file(fastq_files) from TRIMMED_SAMPLES_FOR_FASTQC

  output:
    set file("*_fastqc.html"), file("*_fastqc.zip") into FASTQC_2_OUTPUT
    set val(sample_id), val(1) into CLEAN_TRIMMED_FASTQ_FASTQC_SIGNAL

  script:
  """
  echo "#TRACE sample_id=${sample_id}"
  echo "#TRACE trimmed_fastq_lines=`cat *.fastq | wc -l`"

  fastqc ${fastq_files}
  """
}



/**
 * Performs hisat2 alignment of fastq files to a genome reference
 *
 * depends: trimmomatic
 */
process hisat2 {
  publishDir "${params.outdir}/Samples/${sample_id}", mode: params.publish_dir_mode, pattern: "*.log"
  tag { sample_id }
  label "multithreaded"
  label "hisat2"

  input:
    set val(sample_id), file(fastq_files) from TRIMMED_SAMPLES_FOR_HISAT2
    file indexes from HISAT2_INDEXES

  output:
    set val(sample_id), file("*.sam") into HISAT2_SAM_FILE
    set val(sample_id), file("*.sam.log") into HISAT2_SAM_LOG
    set val(sample_id), file("*.sam") into SAM_FOR_CLEANING
    set val(sample_id), val(1) into CLEAN_TRIMMED_FASTQ_HISAT_SIGNAL
    set val(sample_id), val(1) into CLEAN_MERGED_FASTQ_HISAT_SIGNAL

  script:
  """
  echo "#TRACE sample_id=${sample_id}"
  echo "#TRACE n_cpus=${task.cpus}"
  echo "#TRACE trimmed_fastq_lines=`cat *.fastq | wc -l`"
  echo "#TRACE index_bytes=`stat -Lc '%s' ${indexes} | awk '{sum += \$1} END {print sum}'`"

  hisat2.sh \
    ${sample_id} \
    ${params.hisat2_base_name} \
    ${task.cpus} \
    "${fastq_files}" \
  """
}



/**
 * Sorts the SAM alignment file and coverts it to binary BAM
 *
 * depends: hisat2
 */
process samtools_sort {
  publishDir "${params.outdir}/Samples/${sample_id}", mode: params.publish_dir_mode, pattern: publish_pattern_samtools_sort
  tag { sample_id }
  label "samtools"

  input:
    set val(sample_id), file(sam_file) from HISAT2_SAM_FILE

  output:
    set val(sample_id), file("*.bam") into SORTED_FOR_INDEX
    set val(sample_id), file("*.bam") into BAM_FOR_CLEANING
    set val(sample_id), val(1) into CLEAN_SAM_SIGNAL

  script:
  """
  echo "#TRACE sample_id=${sample_id}"
  echo "#TRACE sam_lines=`cat *.sam | wc -l`"

  samtools sort \
    -o ${sample_id}.bam \
    -O bam \
    -T temp \
    ${sam_file}
  """
}



/**
 * Indexes the BAM alignment file
 *
 * depends: samtools_index
 */
process samtools_index {
  publishDir "${params.outdir}/Samples/${sample_id}", mode: params.publish_dir_mode, pattern: publish_pattern_samtools_index
  tag { sample_id }
  label "samtools"

  input:
    set val(sample_id), file(bam_file) from SORTED_FOR_INDEX

  output:
    set val(sample_id), file(bam_file) into BAM_INDEXED_FOR_STRINGTIE
    set val(sample_id), file("*.bam.bai") into BAI_INDEXED_FILE
    set val(sample_id), file("*.bam.log") into BAM_INDEXED_LOG

  script:
  """
  echo "#TRACE sample_id=${sample_id}"
  echo "#TRACE bam_bytes=`stat -Lc '%s' *.bam`"

  samtools index ${bam_file}
  samtools stats ${bam_file} > ${sample_id}.bam.log
  """
}



/**
 * Generates expression-level transcript abundance
 *
 * depends: samtools_index
 */
process stringtie {
  publishDir "${params.outdir}/Samples/${sample_id}", mode: params.publish_dir_mode, pattern: publish_pattern_stringtie_gtf_and_ga
  tag { sample_id }
  label "multithreaded"
  label "stringtie"

  input:
    set val(sample_id), file(bam_file) from BAM_INDEXED_FOR_STRINGTIE
    file gtf_file from GTF_FILE

  output:
    set val(sample_id), file("*.Hisat2.ga"), file("*.Hisat2.gtf") into STRINGTIE_GTF_FOR_FPKM
    set val(sample_id), file("*.Hisat2.*") into STRINGTIE_GTF_FOR_CLEANING
    set val(sample_id), val(1) into CLEAN_BAM_SIGNAL

  script:
  """
  echo "#TRACE sample_id=${sample_id}"
  echo "#TRACE bam_bytes=`stat -Lc '%s' *.bam`"
  echo "#TRACE gtf_lines=`cat *.gtf | wc -l`"

  stringtie \
    -v \
    -p ${task.cpus} \
    -e \
    -o ${sample_id}.Hisat2.gtf \
    -G ${gtf_file} \
    -A ${sample_id}.Hisat2.ga \
    -l ${sample_id} \
    ${bam_file}
  """
}



/**
 * Generates the final FPKM / TPM / raw files from Hisat2
 */
process hisat2_fpkm_tpm {
  publishDir "${params.outdir}/Samples/${sample_id}", mode: params.publish_dir_mode
  tag { sample_id }
  label "stringtie"

  input:
    set val(sample_id), file(ga_file), file(gtf_file) from STRINGTIE_GTF_FOR_FPKM

  output:
    file "*.Hisat2.fpkm" optional true into HISAT2_FPKM
    file "*.Hisat2.tpm" optional true into HISAT2_TPM
    file "*.Hisat2.raw" optional true into HISAT2_RAW
    set val(sample_id), val(1) into CLEAN_STRINGTIE_SIGNAL
    val sample_id into HISAT2_SAMPLE_COMPLETE_SIGNAL

  script:
  """
  echo "#TRACE sample_id=${sample_id}"
  echo "#TRACE publish_fpkm=${params.hisat2_keep_fpkm}"
  echo "#TRACE publish_tpm=${params.hisat2_keep_tpm}"
  echo "#TRACE publish_raw=${params.hisat2_keep_counts}"
  echo "#TRACE ga_lines=`cat *.ga | wc -l`"
  echo "#TRACE gtf_lines=`cat *.gtf | wc -l`"

  hisat2_fpkm_tpm.sh \
    ${params.hisat2_keep_fpkm} \
    ${sample_id} \
    ${params.hisat2_keep_tpm} \
    ${params.hisat2_keep_counts}
  """
}



/**
 * The multiqc process should run when all samples have
 * completed or if on a resume when the bootstrap signal is
 * received.
 */
MULTIQC_RUN = MULTIQC_READY_SIGNAL.mix(MULTIQC_BOOTSTRAP)

FASTQC_1_OUTPUT.flatten()
  .concat(FASTQC_2_OUTPUT.flatten(), KALLISTO_LOG, SALMON_GA_LOG, HISAT2_SAM_LOG, BAM_INDEXED_LOG, TRIMMED_SAMPLE_LOG).set {MULTIQC_FILES}

/**
 * Process to generate the multiqc report once everything is completed
 */
process multiqc {
  label "multiqc"
  cache false
  publishDir "${params.outdir}/reports", mode: params.publish_dir_mode

  input:
    val signal from MULTIQC_RUN.collect()
    file multiqc_config from MULTIQC_CONFIG
    file gemmaker_logo from MULTIQC_CUSTOM_LOGO
    file input_files from  MULTIQC_FILES.collect()

  output:
    file "multiqc_data" into MULTIQC_DATA
    file "multiqc_report.html" into MULTIQC_REPORT

  when:
    params.publish_multiqc_report == true

  script:
  """
  multiqc --config ${multiqc_config} ./

  """
}



/**
 * The create_gem process should run when all samples have
 * completed or if on a resume when the bootstrap signal is
 * received.
 */
CREATE_GEM_RUN = CREATE_GEM_READY_SIGNAL.mix(CREATE_GEM_BOOTSTRAP)
SALMON_RAW
  .concat(SALMON_TPM, KALLISTO_RAW, KALLISTO_TPM, HISAT2_RAW, HISAT2_TPM, HISAT2_FPKM).set {QUANT_FILES}

/**
 * Creates the GEM file from all the FPKM/TPM outputs
 */
process create_gem {
  label "create_gem"
  publishDir "${params.outdir}/GEMs", mode: params.publish_dir_mode

  input:
    val signal from CREATE_GEM_RUN.collect()
    file input_files from QUANT_FILES.collect()

  output:
    file "*.GEM.*.txt" into GEM_FILES

  when:
    publish_gem == true

  script:
  """
  echo "#TRACE publish_fpkm=${publish_fpkm}"
  echo "#TRACE hisat2_enable=${hisat2_enable}"
  echo "#TRACE publish_tpm=${publish_tpm}"
  echo "#TRACE publish_raw=${publish_raw}"
  echo "#TRACE fpkm_lines=`cat ${params.outdir}/*.fpkm 2> /dev/null  | wc -l`"
  echo "#TRACE tpm_lines=`cat ${params.outdir}/*.tpm 2> /dev/null | wc -l`"
  echo "#TRACE raw_lines=`cat ${params.outdir}/*.raw 2> /dev/null | wc -l`"

  create_gem.sh \
    ${publish_fpkm} \
    ${hisat2_enable} \
    . \
    GEMmaker \
    ${publish_raw} \
    ${publish_tpm}
  """
}

/**
 * Creates a report of any SRA run IDs that failed and why they failed.
 */
process failed_run_report {
  label "reports"
  publishDir "${params.outdir}/reports", mode: params.publish_dir_mode

  input:
    file metadata_failed_runs from METADATA_FAILED_RUNS
    file download_failed_runs from DOWNLOAD_FAILED_RUNS.collect()
    file fastq_dump_failed_runs from FASTQ_DUMP_FAILED_RUNS.collect()
    file failed_run_template from FAILED_RUN_TEMPLATE

  output:
    file "failed_SRA_run_report.html" into FAILED_RUN_REPORT

  script:
  """
  failed_runs_report.py --template ${failed_run_template}
  """

}



/**
 * PROCESSES FOR CLEANING LARGE FILES
 *
 * Nextflow doesn't allow files to be removed from the
 * work directories that are used in Channels.  If it
 * detects a different timestamp or change in file
 * size than what was cached it will rerun the process.
 * To trick Nextflow we will truncate the file to a
 * sparce file of size zero but masquerading as its
 * original size, we will also reset the original modify
 * and access times.
 */

/**
 * Cleans downloaded SRA files
 */
SRA_TO_CLEAN
  .mix(CLEAN_SRA_SIGNAL)
  .groupTuple(size: 2)
  .set { CLEAN_SRA_READY }

process clean_sra {
  tag { sample_id }

  input:
    set val(sample_id), val(files_list) from CLEAN_SRA_READY

  when:
    params.keep_sra == false

  script:
  """
  clean_work_files.sh "${files_list[0]}"
  """
}

FCLEAN = DOWNLOADED_FASTQ_FOR_CLEANING.mix(CLEAN_FASTQ_SIGNAL)
FCLEAN.groupTuple(size: 2).set { FASTQ_CLEANUP_READY }

process clean_fastq {
  tag { sample_id }

  input:
    set val(sample_id), val(files_list) from FASTQ_CLEANUP_READY

  when:
    params.keep_retrieved_fastq == false

  script:
  flist = files_list[0].join(" ")
  """
  clean_work_files.sh "${flist}"
  """
}

/**
 * The "signal" to clean up the merged FASTQ files consists of:
 *
 *   1. The list of FASTQ files in the work directory to clean up. This comes
 *      from the MERGED_FASTQ_FOR_CLEANING channel.
 *   2. The signal from the alignment tool that it is complete. This comes
 *      from one of these channels:  CLEAN_MERGED_FASTQ_HISAT_SIGNAL,
 *      CLEAN_MERGED_FASTQ_KALLISTO_SIGNAL, CLEAN_MERGED_FASTQ_SALMON_SIGNAL.
 *   3. The signal from the FastQC tool that it is complete. This comes
 *      from this signal: CLEAN_MERGED_FASTQ_FASTQC_SIGNAL.
 *
 * Each element from these channels is a list (or tuple), where the key
 * (first element) is the sample ID.  We use a groupTuple to ensure that
 * three elements are combined together into a single tuple that gets emited
 * into the MERGED_FASTQ_CLEANUP_READY channel once all 3 signals are received.
 * Only when all the tuple has all three values is it emited.
 *
 * The first element of the tuple should be the files because they should
 * always come first via the MERGED_FASTQ_FOR_CLEANING channel.
 */
MFCLEAN = MERGED_FASTQ_FOR_CLEANING.mix(CLEAN_MERGED_FASTQ_HISAT_SIGNAL, CLEAN_MERGED_FASTQ_KALLISTO_SIGNAL, CLEAN_MERGED_FASTQ_SALMON_SIGNAL, CLEAN_MERGED_FASTQ_FASTQC_SIGNAL)
MFCLEAN.groupTuple(size: 3).set { MERGED_FASTQ_CLEANUP_READY }

/**
 * Cleans merged fastq files
 */
process clean_merged_fastq {
  tag { sample_id }

  input:
    set val(sample_id), val(files_list) from MERGED_FASTQ_CLEANUP_READY

  when:
    params.keep_retrieved_fastq == false

  script:
  flist = files_list[0].join(" ")
  """
  clean_work_files.sh "${flist}"
  """
}

/**
 * Merge the Trimmomatic samples with Hisat's signal and FastQC signal. Once
 * both tools send the signal that they are done with the trimmed file it can
 * be removed.
 */
TRHIMIX = TRIMMED_FASTQ_FOR_CLEANING.mix(CLEAN_TRIMMED_FASTQ_HISAT_SIGNAL, CLEAN_TRIMMED_FASTQ_FASTQC_SIGNAL)
TRHIMIX.groupTuple(size: 3).set { TRIMMED_FASTQ_CLEANUP_READY }

/**
 * Cleans trimmed fastq files
 */
process clean_trimmed_fastq {
  tag { sample_id }

  input:
    set val(sample_id), val(files_list) from TRIMMED_FASTQ_CLEANUP_READY

  when:
    params.trimmomatic_keep_trimmed_fastq == false

  script:
  flist = files_list[0].join(" ")
  """
  clean_work_files.sh "${flist}"
  """
}

/**
 * Merge the HISAT sam file with samtools_sort signal that it is
 * done so that we can remove these files.
 */
HISSMIX = SAM_FOR_CLEANING.mix(CLEAN_SAM_SIGNAL)
HISSMIX.groupTuple(size: 2).set { SAM_CLEANUP_READY }

/**
 * Clean up SAM files
 */
process clean_sam {
  tag { sample_id }

  input:
    set val(sample_id), val(files_list) from SAM_CLEANUP_READY

  when:
    params.hisat2_keep_sam == false

  script:
  """
  clean_work_files.sh "${files_list[0]}"
  """
}

/**
 * Merge the samtools_sort bam file with stringtie signal that it is
 * done so that we can remove these files.
 */
SSSTMIX = BAM_FOR_CLEANING.mix(CLEAN_BAM_SIGNAL)
SSSTMIX.groupTuple(size: 2).set { BAM_CLEANUP_READY }

/**
 * Clean up BAM files
 */
process clean_bam {
  tag { sample_id }

  input:
    set val(sample_id), val(files_list) from BAM_CLEANUP_READY

  when:
    params.hisat2_keep_bam == false

  script:
  """
  clean_work_files.sh "${files_list[0]}"
  """
}

/**
 * Merge the Kallisto .ga file with the clean kallisto ga signal so that we can
 * remove the .ga file after it has been used
 */
KGAMIX = KALLISTO_GA_TO_CLEAN.mix(CLEAN_KALLISTO_GA_SIGNAL)
KGAMIX.groupTuple(size: 2).set { KALLISTO_GA_CLEANUP_READY }

/**
 * Clean up Kallisto GA files
 */
process clean_kallisto_ga {
  tag { sample_id }

  input:
    set val(sample_id), val(directory) from KALLISTO_GA_CLEANUP_READY

  when:
    params.kallisto_keep_data == false

  script:
  """
  clean_work_dirs.sh "${directory[0]}"
  """
}

/**
 * Merge the Salmon .ga file with the clean salmon ga signal so that we can
 * remove the .ga file after it has been used
 */
SGAMIX = SALMON_GA_TO_CLEAN.mix(CLEAN_SALMON_GA_SIGNAL)
SGAMIX.groupTuple(size: 2).set { SALMON_GA_CLEANUP_READY }

/**
 * Clean up Salmon GA files
 */
process clean_salmon_ga {
  tag { sample_id }

  input:
    set val(sample_id), val(files_list) from SALMON_GA_CLEANUP_READY

  when:
    params.salmon_keep_data == false

  script:
  """
  clean_work_files.sh "${files_list[0]}"
  """
}

/**
 * Merge the Salmon .ga file with the clean salmon ga signal so that we can
 * remove the .ga file after it has been used
 */
SMIX = STRINGTIE_GTF_FOR_CLEANING.mix(CLEAN_STRINGTIE_SIGNAL)
SMIX.groupTuple(size: 2).set { STRINGTIE_CLEANUP_READY }

/**
 * Clean up Salmon GA files
 */
process clean_stringtie_ga {
  tag { sample_id }

  input:
    set val(sample_id), val(files_list) from STRINGTIE_CLEANUP_READY

  when:
    params.hisat2_keep_data == false

  script:
  flist = files_list[0].join(" ")
  """
  clean_work_files.sh "${flist}"
  """
}
