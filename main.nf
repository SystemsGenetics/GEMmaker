#!/usr/bin/env nextflow

import java.nio.channels.FileLock
import java.nio.channels.FileChannel
import java.nio.channels.OverlappingFileLockException;

/**
 * ========
 * GEMmaker
 * ========
 *
 * Authors:
 *  + John Hadish
 *  + Tyler Biggs
 *  + Stephen Ficklin
 *  + Ben Shealy
 *  + Connor Wytko
 *
 * Summary:
 *   A workflow for processing a large amount of RNA-seq data
 */

println """\

===================================
 G E M M A K E R   P I P E L I N E
===================================

Workflow Information:
---------------------
  Project Directory:  ${workflow.projectDir}
  Launch Directory:   ${workflow.launchDir}
  Work Directory:     ${workflow.workDir}
  Config Files:       ${workflow.configFiles}
  Container Engine:   ${workflow.containerEngine}
  Profile(s):         ${workflow.profile}


Input Parameters:
-----------------
  Remote fastq list path:     ${params.input.input_data_dir}/${params.input.remote_sample_list}
  Local sample glob:          ${params.input.input_data_dir}/${params.input.local_sample_files}


Quantification Tool Input:
--------------------------
  Use Hisat2:                 ${params.input.hisat2.enable}
  Hisat2 Index Prefix:        ${params.input.reference_name}
  Hisat2 GTF File:            ${params.input.hisat2.gtf_file}
  Hisat2 Index Directory:     ${params.input.hisat2.index_dir}

  Use Kallisto:               ${params.input.kallisto.enable}
  Kallisto Index File:        ${params.input.kallisto.index_file}

  Use Salmon:                 ${params.input.salmon.enable}
  Salmon Index Directory:     ${params.input.salmon.index_dir}


Output Parameters:
------------------
  Output directory:           ${params.output.dir}
  Publish SRA:                ${params.output.publish_sra}
  Publish downloaded FASTQ:   ${params.output.publish_downloaded_fastq}
  Publish trimmed FASTQ:      ${params.output.publish_trimmed_fastq}
  Publish BAM:                ${params.output.publish_bam}
  Publish Gene Abundance:     ${params.output.publish_gene_abundance}
  Publish GTF_GA:             ${params.output.publish_stringtie_gtf_and_ga}
  Publish RAW:                ${params.output.publish_raw}
  Publish FPKM:               ${params.output.publish_fpkm}
  Publish TPM:                ${params.output.publish_tpm}
  MultiQC:                    ${params.output.multiqc}
  Create GEM:                 ${params.output.create_gem}


Execution Parameters:
---------------------
  Queue size:                 ${params.execution.queue_size}


Software Parameters:
--------------------
  Trimmomatic clip path:      ${params.software.trimmomatic.clip_path}
  Trimmomatic minimum ratio:  ${params.software.trimmomatic.MINLEN}
"""



// Indicates if a tool was selected.
has_tool = 0

// Indicates which tool the user selected.
selected_tool = 0

// Print out details per the selected tool.
if (params.input.hisat2.enable == true) {
  has_tool++
  selected_tool = 0
}
if (params.input.kallisto.enable == true) {
  has_tool++
  selected_tool = 1
}
if (params.input.salmon.enable == true) {
  has_tool++
  selected_tool = 2
}

if (has_tool == 0) {
  error "Error: You must select a valid quantification tool in the 'nextflow.config' file"
}
if (has_tool > 1) {
  error "Error: Please select only one quantification tool in the 'nextflow.config' file"
}

// Check to make sure that required reference files exist
// If Hisat2 was selected:
if (selected_tool == 0)
{
  gtfFile = file("${params.input.reference_dir}/${params.input.hisat2.gtf_file}")
  if (gtfFile.isEmpty())
  {
    error "Error: GTF reference file for Hisat2 does not exist or is empty! Please Check that you have the proper references, that they are placed in the reference directory, and they are named properly.\
    \nGEMmaker is missing the following file: '${params.input.reference_dir}/${params.input.hisat2.gtf_file}' (where '*' is the name of your organism)"
  }
  hisat2_index_dir = file("${params.input.reference_dir}/${params.input.hisat2.index_dir}")
  if(!hisat2_index_dir.isDirectory())
  {
    error "Error: hisat2 Index Directory does not exist or is empty! Please Check that you have the proper references, that they are placed in the reference directory, and they are named properly.\
    \nGEMmaker is missing the following file: '${params.input.reference_dir}/${params.input.hisat2.index_dir}' (where '*' is the name of your organism)"
  }

}
// If Kallisto was selected
if (selected_tool == 1)
{
  kallisto_index_file = file("${params.input.reference_dir}/${params.input.kallisto.index_file}")
  if (kallisto_index_file.isEmpty())
  {
    error "Error: Kallisto Index File does not exist or is empty! Please Check that you have the proper references, that they are placed in the reference directory, and they are named properly.\
    \nGEMmaker is missing the following file: '${params.input.reference_dir}/${params.input.kallisto.index_file}' (where '*' is the name of your organism)"
  }
}
// If Salmon was selected
if (selected_tool == 2)
{
  salmon_index_dir = file("${params.input.reference_dir}/${params.input.salmon.index_dir}")
  if (!salmon_index_dir.isDirectory())
  {
    error "Error: Salmon Index Directory does not exist or is empty! Please Check that you have the proper references, that they are placed in the reference directory, and they are named properly.\
    \nGEMmaker is missing the following file: '${params.input.reference_dir}/${params.input.salmon.index_dir}' (where '*' is the name of your organism)"
  }
}

/**
 * Create value channels that can be reused
 */
HISAT2_INDEXES = Channel.fromPath("${params.input.reference_dir}/${params.input.hisat2.index_dir}/*").collect()
KALLISTO_INDEX = Channel.fromPath("${params.input.reference_dir}/${params.input.kallisto.index_file}").collect()
SALMON_INDEXES = Channel.fromPath("${params.input.reference_dir}/${params.input.salmon.index_dir}/*").collect()
FASTA_ADAPTER = Channel.fromPath("${params.software.trimmomatic.clip_path}").collect()
GTF_FILE = Channel.fromPath("${params.input.reference_dir}/${params.input.hisat2.gtf_file}").collect()



/**
 * Local Sample Input.FASTA_ADAPTER
 * This checks the folder that the user has given
 */
if (params.input.local_sample_files == "none") {
  Channel.empty().set { LOCAL_SAMPLE_FILES_FOR_STAGING }
  Channel.empty().set { LOCAL_SAMPLE_FILES_FOR_JOIN }
}
else {
  Channel.fromFilePairs( "${params.input.input_data_dir}/${params.input.local_sample_files}", size: -1 )
    .set { LOCAL_SAMPLE_FILES_FOR_STAGING }
  Channel.fromFilePairs( "${params.input.input_data_dir}/${params.input.local_sample_files}", size: -1 )
    .set { LOCAL_SAMPLE_FILES_FOR_JOIN }
}

/**
 * Remote fastq_run_id Input.
 */
if (params.input.remote_sample_list == "none") {
  Channel.empty().set { SRR_FILE }
}
else {
  Channel.fromPath("${params.input.input_data_dir}/${params.input.remote_sample_list}").set { SRR_FILE }
}

/**
 * Make sure that at least one output format is enabled.
 */
if ( params.input.hisat2.enabled == true && params.output.publish_raw == false && params.output.publish_fpkm == false && params.output.publish_tpm == false ) {
  error "Error: at least one output format (raw, fpkm, tpm) must be enabled for hisat2"
}

if ( params.input.hisat2.enabled == false && params.output.publish_raw == false && params.output.publish_tpm == false ) {
  error "Error: at least one output format (raw, tpm) must be enabled for kallisto / salmon"
}

/**
 * Set the pattern for publishing downloaded FASTQ files
 */
publish_pattern_fastq_dump = "{none}";
if (params.output.publish_downloaded_fastq == true) {
  publish_pattern_fastq_dump = "{*.fastq}";
}



/**
 * Set the pattern for publishing trimmed FASTQ files
 */
publish_pattern_trimmomatic = "{*.trim.log}";
if (params.output.publish_trimmed_fastq == true) {
  publish_pattern_trimmomatic = "{*.trim.log,*_trim.fastq}";
}



/**
 * Set the pattern for publishing BAM files
 */
publish_pattern_samtools_sort = "{*.log}";
publish_pattern_samtools_index = "{*.log}";

if (params.output.publish_bam == true) {
  publish_pattern_samtools_sort = "{*.log,*.bam}";
  publish_pattern_samtools_index = "{*.log,*.bam.bai}";
}

/**
 * Set the pattern for publishing Kallisto GA files
 */
publish_pattern_Kallisto_GA = "{*.log}";
if (params.output.publish_gene_abundance == true) {
  publish_pattern_Kallisto_GA = "{*.ga,*.log}";
}

/**
 * Set the pattern for publishing Salmon GA files
 * Publishs only log file used by multiqc if false
 */
publish_pattern_Salmon_GA = "{*.ga/aux_info/meta_info.json,*.ga/libParams/flenDist.txt}"
if (params.output.publish_gene_abundance == true) {
  publish_pattern_Salmon_GA = "{*.ga}";
}


/**
 * Set the pattern for publishing STRINGTIE GA and GTF files
 */
publish_pattern_stringtie_gtf_and_ga = "{none}"
if (params.output.publish_stringtie_gtf_and_ga == true) {
  publish_pattern_stringtie_gtf_and_ga = "{*.ga, *.gtf}";
}


/**
 * Make sure that at least one output format is enabled.
 */
if ( params.software.alignment == 0 && params.output.publish_raw == false && params.output.publish_fpkm == false && params.output.publish_tpm == false ) {
  error "Error: at least one output format (raw, fpkm, tpm) must be enabled for hisat2"
}

if ( params.software.alignment != 0 && params.output.publish_raw == false && params.output.publish_tpm == false ) {
  error "Error: at least one output format (raw, tpm) must be enabled for kallisto / salmon"
}



/**
 * Retrieves metadata for all of the remote samples
 * and maps SRA runs to SRA experiments.
 */
process retrieve_sra_metadata {
  publishDir params.output.dir, mode: params.output.publish_mode, pattern: "*.GEMmaker.meta.*", saveAs: { "${it.tokenize(".")[0]}/${it}" }
  label "python3"

  input:
    file srr_file from SRR_FILE

  output:
    stdout REMOTE_SAMPLES_LIST
    file "*.GEMmaker.meta.*"

  script:
    """
    retrieve_sra_metadata.py ${srr_file}
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
  .map{ [it[0], it[1], 'local' ] }
  .set{LOCAL_SAMPLES_FOR_STAGING}

ALL_SAMPLES = REMOTE_SAMPLES_FOR_STAGING
  .mix(LOCAL_SAMPLES_FOR_STAGING)

// Create the directories we'll use for running
// batches
file("${workflow.workDir}/GEMmaker").mkdir()
file("${workflow.workDir}/GEMmaker/stage").mkdir()
file("${workflow.workDir}/GEMmaker/process").mkdir()
file("${workflow.workDir}/GEMmaker/done").mkdir()

// Channels to bootstrap post-processing of
// sample results if a resume is performed when
// all samples have completed.
MULTIQC_BOOTSTRAP = Channel.create()
CREATE_GEM_BOOTSTRAP = Channel.create()

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
  executor "local"
  tag {sample[0]}

  input:
    val sample from ALL_SAMPLES

  output:
    val (1) into SAMPLES_READY_SIGNAL

  exec:
    // Get any samples to skip
    skip_samples = []
    skip_file = file("${params.input.input_data_dir}/${params.input.skip_list_path}")
    if (skip_file.exists()) {
      skip_file.eachLine { line ->
        skip_samples << line.trim()
      }
    }

    // Only stage files that should not be skipped.
    if (skip_samples.intersect([sample[0]]) == [])  {

      // Create a file for each samples.
      sample_file = file("${workflow.workDir}/GEMmaker/stage/" + sample[0] + '.sample.csv')
      sample_file.withWriter {

        // Get the sample type: local or remote.
        type = sample[2]

        // If this is a local file.
        if (type.equals('local')) {
          if (sample[1].size() > 1) {
            files = sample[1]
            files_str = files.join('::')
            it.writeLine '"' + sample[0] + '","' + files_str + '","' + type + '"'
          }
          else {
            it.writeLine '"' + sample[0] + '","' + sample[1].first().toString() + '","' + type + '"'
          }
        }
        // If this is a remote file.
        else {
          it.writeLine '"' + sample[0] + '","' + sample[1] + '","' + type + '"'
        }
      }
    }
}



// When all batch files are created we need to then
// move the first file into the process directory.
SAMPLES_READY_SIGNAL.collect().set { FIRST_SAMPLE_START_SIGNAL }



/**
 * Moves the first set of sample files into the process directory.
 */
process start_first_batch {
  executor "local"
  cache false

  input:
    val signal from FIRST_SAMPLE_START_SIGNAL

  exec:
    // Move the first set of sample file into the processing directory
    // so that we jumpstart the workflow.
    sample_files = file("${workflow.workDir}/GEMmaker/stage/*.sample.csv");
    start_samples = sample_files.sort().take(params.execution.queue_size)
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
  executor "local"
  cache false
  tag { sample_file }

  input:
    file(sample_file) from NEXT_SAMPLE

  output:
    stdout SAMPLE_FILE_CONTENTS

  script:
    """
      cat $sample_file
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

// Create the channel that will collate all the signals
// and release a signal when the sample is complete
SAMPLE_COMPLETE_SIGNAL = Channel.create()
SAMPLE_COMPLETE_SIGNAL
  .mix(HISAT2_SAMPLE_COMPLETE_SIGNAL, KALLISTO_SAMPLE_COMPLETE_SIGNAL, SALMON_SAMPLE_COMPLETE_SIGNAL)
  .into { NEXT_SAMPLE_SIGNAL; MULTIQC_READY_SIGNAL; CREATE_GEM_READY_SIGNAL }



/**
 * Handles the end of a sample by moving a new sample
 * file into the process directory which triggers
 * the NEXT_SAMPLE.watchPath channel.
 */
process next_sample {
  executor "local"

  input:
    val sample_id from NEXT_SAMPLE_SIGNAL

  exec:

    // Use a file lock to prevent a race condition for grabbing the next sample.
    File lockfile = null;
    FileChannel channel = null;
    FileLock lock = null;
    success = false

    try {
      attempts = 0
      while (!lock)  {
        // We will try for a release lock for about 10 minutes max. The sleep
        // time is set to 1 second.  It may take this long on a slow NFS
        // system.
        if (attempts < 6000) {
          try {
            lockfile = new File("${workflow.workDir}/GEMmaker/gemmaker.lock")
            channel = new RandomAccessFile(lockfile, "rw").getChannel()
            lock = channel.lock()
          }
          catch (OverlappingFileLockException e) {
            // Do nothing, let's try a few more times....
          }
          if (!lock) {
            println "Waiting on lock. attempt " + attempts + "..."
            sleep 1000
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
        lock.release();
      }
      if (channel) {
        channel.close();
      }
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
  tag { sample_id }
  label "sratoolkit"

  input:
    set val(sample_id), val(run_ids), val(type) from REMOTE_SAMPLES

  output:
    set val(sample_id), file("*.sra") into SRA_TO_EXTRACT
    set val(sample_id), file("*.sra") into SRA_TO_CLEAN

  script:
  """
    ids=`echo $run_ids | perl -p -e 's/[\\[,\\]]//g' | perl -p -e 's/\\s*\$//' | perl -p -e 's/\\s+/,/g'`
    retrieve_sra.py \$ids
  """
}



/**
 * Extracts FASTQ files from downloaded SRA files.
 */
process fastq_dump {
  publishDir params.output.dir, mode: params.output.publish_mode, pattern: publish_pattern_fastq_dump, saveAs: { "${sample_id}/${it}" }
  tag { sample_id }
  label "sratoolkit"

  input:
    set val(sample_id), val(sra_files) from SRA_TO_EXTRACT

  output:
    set val(sample_id), file("*.fastq") into DOWNLOADED_FASTQ_FOR_MERGING
    set val(sample_id), file("*.fastq") into DOWNLOADED_FASTQ_FOR_CLEANING
    set val(sample_id), val(1) into CLEAN_SRA_SIGNAL

  script:
  """
  files=`echo $sra_files | perl -p -e 's/[\\[,\\]]//g' | perl -p -e 's/\\s*\$//' | perl -p -e 's/\\s+/,/g'`
  sra2fastq.py \$files
  """
}



/**
 * This process merges the fastq files based on their sample_id number.
 */
process fastq_merge {
  tag { sample_id }

  input:
    set val(sample_id), file(grouped) from DOWNLOADED_FASTQ_FOR_MERGING

  output:
    set val(sample_id), file("${sample_id}_?.fastq") into MERGED_SAMPLES_FOR_COUNTING
    set val(sample_id), file("${sample_id}_?.fastq") into MERGED_SAMPLES_FOR_FASTQC_1
    set val(sample_id), file("${sample_id}_?.fastq") into MERGED_FASTQ_FOR_CLEANING

  /**
   * This command tests to see if ls produces a 0 or not by checking
   * its standard out. We do not use a "if [-e *foo]" becuase it gets
   * confused if there are more than one things returned by the wildcard
   */
  script:
  """
  fastq_merge.sh ${sample_id} ${params.input.publish_downloaded_fastq}
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
  publishDir params.output.sample_dir, mode: params.output.publish_mode, pattern: "*_fastqc.*"
  tag { sample_id }
  label "fastqc"

  input:
    set val(sample_id), file(pass_files) from COMBINED_SAMPLES_FOR_FASTQC_1

  output:
    set file("${sample_id}_?_fastqc.html") , file("${sample_id}_?_fastqc.zip") optional true into FASTQC_1_OUTPUT
    set val(sample_id), val(1) into CLEAN_MERGED_FASTQ_FASTQC_SIGNAL

  script:
  """
  fastqc $pass_files
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
  publishDir params.output.sample_dir, mode: params.output.publish_mode, pattern: publish_pattern_Kallisto_GA
  tag { sample_id }
  label "kallisto"

  input:
    set val(sample_id), file(pass_files) from KALLISTO_CHANNEL
    file kallisto_index from KALLISTO_INDEX

  output:
    set val(sample_id), file("${sample_id}_vs_${params.input.reference_name}.Kallisto.ga") into KALLISTO_GA
    set val(sample_id), file("${sample_id}_vs_${params.input.reference_name}.Kallisto.ga") into KALLISTO_GA_TO_CLEAN
    set val(sample_id), val(1) into CLEAN_MERGED_FASTQ_KALLISTO_SIGNAL
    file "*kallisto.log" into KALLISTO_LOG

  script:
  """
  kallisto.sh ${sample_id} ${kallisto_index} ${params.input.reference_name}
  """
}



/**
 * Generates the final TPM and raw count files for Kallisto
 */
process kallisto_tpm {
  publishDir params.output.sample_dir, mode: params.output.publish_mode
  tag { sample_id }

  input:
    set val(sample_id), file("${sample_id}_vs_${params.input.reference_name}.Kallisto.ga") from KALLISTO_GA

  output:
    file "${sample_id}_vs_${params.input.reference_name}.Kallisto.tpm" optional true into KALLISTO_TPM
    file "${sample_id}_vs_${params.input.reference_name}.Kallisto.raw" optional true into KALLISTO_RAW
    set val(sample_id), val(1) into CLEAN_KALLISTO_GA_SIGNAL
    val sample_id  into KALLISTO_SAMPLE_COMPLETE_SIGNAL

  script:
  """
  kallisto_tpm.sh ${sample_id} ${params.output.publish_tpm} ${params.output.publish_raw} ${params.input.reference_name}
  """
}



/**
 * Performs SALMON alignemnt of fastq files
 */
process salmon {
  publishDir params.output.sample_dir, mode: params.output.publish_mode, pattern: publish_pattern_Salmon_GA
  tag { sample_id }
  label "multithreaded"
  label "salmon"

  input:
    set val(sample_id), file(pass_files) from SALMON_CHANNEL
    file salmon_index from SALMON_INDEXES

  output:
    set val(sample_id), file("${sample_id}_vs_${params.input.reference_name}.Salmon.ga") into SALMON_GA
    set val(sample_id), file("*.ga/aux_info/meta_info.json"), file("*.ga/libParams/flenDist.txt") into SALMON_GA_LOG
    set val(sample_id), file("${sample_id}_vs_${params.input.reference_name}.Salmon.ga/quant.sf") into SALMON_GA_TO_CLEAN
    set val(sample_id), val(1) into CLEAN_MERGED_FASTQ_SALMON_SIGNAL

  script:
  """
  salmon.sh ${sample_id} ${task.cpus} ${params.input.reference_name}
  """
}



/**
 * Generates the final TPM file for Salmon
 */
process salmon_tpm {
  publishDir params.output.sample_dir, mode: params.output.publish_mode
  tag { sample_id }

  input:
    set val(sample_id), file("${sample_id}_vs_${params.input.reference_name}.Salmon.ga") from SALMON_GA

  output:
    file "${sample_id}_vs_${params.input.reference_name}.Salmon.tpm" optional true into SALMON_TPM
    file "${sample_id}_vs_${params.input.reference_name}.Salmon.raw" optional true into SALMON_RAW
    set val(sample_id), val(1) into CLEAN_SALMON_GA_SIGNAL
    val sample_id  into SALMON_SAMPLE_COMPLETE_SIGNAL

  script:
  """
  salmon_tpm.sh ${params.output.publish_tpm} ${sample_id} \
  ${params.input.reference_name} ${params.output.publish_raw}
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
 * read length. The percenage is determined by the user in the
 * "nextflow.config" file
 */
process trimmomatic {
  publishDir params.output.sample_dir, mode: params.output.publish_mode, pattern: publish_pattern_trimmomatic
  tag { sample_id }
  label "multithreaded"
  label "trimmomatic"

  input:
    set val(sample_id), file("${sample_id}_?.fastq") from HISAT2_CHANNEL
    file fasta_adapter from FASTA_ADAPTER

  output:
    set val(sample_id), file("${sample_id}_*trim.fastq") into TRIMMED_SAMPLES_FOR_FASTQC
    set val(sample_id), file("${sample_id}_*trim.fastq") into TRIMMED_SAMPLES_FOR_HISAT2
    set val(sample_id), file("${sample_id}_*trim.fastq") into TRIMMED_FASTQ_FOR_CLEANING
    set val(sample_id), file("${sample_id}.trim.log") into TRIMMED_SAMPLE_LOG

  script:
  """
  trimmomatic.sh ${sample_id} ${params.software.trimmomatic.MINLEN} ${task.cpus} ${fasta_adapter} ${params.software.trimmomatic.LEADING} ${params.software.trimmomatic.TRAILING} ${params.software.trimmomatic.SLIDINGWINDOW} ${params.software.trimmomatic.quality}
  """
}

/**
 * Performs fastqc on fastq files post trimmomatic
 * Files are stored to an independent folder
 */
process fastqc_2 {
  publishDir params.output.sample_dir, mode: params.output.publish_mode, pattern: "*_fastqc.*"
  tag { sample_id }
  label "fastqc"

  input:
    set val(sample_id), file(pass_files) from TRIMMED_SAMPLES_FOR_FASTQC

  output:
    set file("${sample_id}_??_trim_fastqc.html"), file("${sample_id}_??_trim_fastqc.zip") optional true into FASTQC_2_OUTPUT
    set val(sample_id), val(1) into CLEAN_TRIMMED_FASTQ_FASTQC_SIGNAL

  script:
  """
  fastqc $pass_files
  """
}



/**
 * Performs hisat2 alignment of fastq files to a genome reference
 *
 * depends: trimmomatic
 */
process hisat2 {
  publishDir params.output.sample_dir, mode: params.output.publish_mode, pattern: "*.log"
  tag { sample_id }
  label "multithreaded"
  label "hisat2"

  input:
    set val(sample_id), file(input_files) from TRIMMED_SAMPLES_FOR_HISAT2
    file indexes from HISAT2_INDEXES

  output:
    set val(sample_id), file("${sample_id}_vs_${params.input.reference_name}.sam") into INDEXED_SAMPLES
    set val(sample_id), file("${sample_id}_vs_${params.input.reference_name}.sam.log") into INDEXED_SAMPLES_LOG
    set val(sample_id), file("${sample_id}_vs_${params.input.reference_name}.sam") into SAM_FOR_CLEANING
    set val(sample_id), val(1) into CLEAN_TRIMMED_FASTQ_HISAT_SIGNAL
    set val(sample_id), val(1) into CLEAN_MERGED_FASTQ_HISAT_SIGNAL

  script:
  """
  hisat2.sh ${sample_id} ${params.input.reference_name} ${task.cpus}
  """
}



/**
 * Sorts the SAM alignment file and coverts it to binary BAM
 *
 * depends: hisat2
 */
process samtools_sort {
  publishDir params.output.sample_dir, mode: params.output.publish_mode, pattern: publish_pattern_samtools_sort
  tag { sample_id }
  label "samtools"

  input:
    set val(sample_id), file("${sample_id}_vs_${params.input.reference_name}.sam") from INDEXED_SAMPLES

  output:
    set val(sample_id), file("${sample_id}_vs_${params.input.reference_name}.bam") into SORTED_FOR_INDEX
    set val(sample_id), file("${sample_id}_vs_${params.input.reference_name}.bam") into BAM_FOR_CLEANING
    set val(sample_id), val(1) into CLEAN_SAM_SIGNAL

  script:
    """
    samtools sort \
      -o ${sample_id}_vs_${params.input.reference_name}.bam \
      -O bam \
      -T temp \
      ${sample_id}_vs_${params.input.reference_name}.sam
    """
}



/**
 * Indexes the BAM alignment file
 *
 * depends: samtools_index
 */
process samtools_index {
  publishDir params.output.sample_dir, mode: params.output.publish_mode, pattern: publish_pattern_samtools_index
  tag { sample_id }
  label "samtools"

  input:
    set val(sample_id), file("${sample_id}_vs_${params.input.reference_name}.bam") from SORTED_FOR_INDEX

  output:
    set val(sample_id), file("${sample_id}_vs_${params.input.reference_name}.bam") into BAM_INDEXED_FOR_STRINGTIE
    set val(sample_id), file("${sample_id}_vs_${params.input.reference_name}.bam.bai") into BAI_INDEXED_FILE
    set val(sample_id), file("${sample_id}_vs_${params.input.reference_name}.bam.log") into BAM_INDEXED_LOG

  script:
    """
    samtools index ${sample_id}_vs_${params.input.reference_name}.bam
    samtools stats ${sample_id}_vs_${params.input.reference_name}.bam > ${sample_id}_vs_${params.input.reference_name}.bam.log
    """
}



/**
 * Generates expression-level transcript abundance
 *
 * depends: samtools_index
 */
process stringtie {
  publishDir params.output.sample_dir, mode: params.output.publish_mode, pattern: publish_pattern_stringtie_gtf_and_ga
  tag { sample_id }
  label "multithreaded"
  label "stringtie"

  input:
    // We don't really need the .bam file, but we want to ensure
    // this process runs after the samtools_index step so we
    // require it as an input file.
    set val(sample_id), file("${sample_id}_vs_${params.input.reference_name}.bam") from BAM_INDEXED_FOR_STRINGTIE
    file gtf_file from GTF_FILE

  output:
    set val(sample_id), file("${sample_id}_vs_${params.input.reference_name}.Hisat2.ga"), file("${sample_id}_vs_${params.input.reference_name}.Hisat2.gtf") into STRINGTIE_GTF_FOR_FPKM
    set val(sample_id), file("${sample_id}_vs_${params.input.reference_name}.Hisat2.*") into STRINGTIE_GTF_FOR_CLEANING
    set val(sample_id), val(1) into CLEAN_BAM_SIGNAL

  script:
    """
    stringtie \
      -v \
      -p ${task.cpus} \
      -e \
      -o ${sample_id}_vs_${params.input.reference_name}.Hisat2.gtf \
      -G ${gtf_file} \
      -A ${sample_id}_vs_${params.input.reference_name}.Hisat2.ga \
      -l ${sample_id} ${sample_id}_vs_${params.input.reference_name}.bam
    """
}


/**
 * Generates the final FPKM / TPM / raw files from Hisat2
 */
process hisat2_fpkm_tpm {
  publishDir params.output.sample_dir, mode: params.output.publish_mode
  tag { sample_id }
  label "stringtie"

  input:
  set val(sample_id), file("${sample_id}_vs_${params.input.reference_name}.Hisat2.ga"), file("${sample_id}_vs_${params.input.reference_name}.Hisat2.gtf") from STRINGTIE_GTF_FOR_FPKM


  output:
    file "${sample_id}_vs_${params.input.reference_name}.Hisat2.fpkm" optional true into FPKMS
    file "${sample_id}_vs_${params.input.reference_name}.Hisat2.tpm" optional true into TPM
    file "${sample_id}_vs_${params.input.reference_name}.Hisat2.raw" optional true into RAW_COUNTS
    set val(sample_id), val(1) into CLEAN_STRINGTIE_SIGNAL
    val sample_id into HISAT2_SAMPLE_COMPLETE_SIGNAL

  script:
  """
  hisat2_fpkm_tpm.sh ${params.output.publish_fpkm} ${sample_id} ${params.input.reference_name} ${params.output.publish_tpm} ${params.output.publish_raw}
  """
}



/**
 * The multiqc process should run when all samples have
 * completed or if on a resume when the bootstrap signal is
 * received.
 */
MULTIQC_RUN = MULTIQC_READY_SIGNAL.mix(MULTIQC_BOOTSTRAP)

/**
 * Process to generate the multiqc report once everything is completed
 */
process multiqc {
  label "multiqc"
  publishDir "${params.output.dir}/reports", mode: params.output.publish_mode

  input:
    val signal from MULTIQC_RUN.collect()

  output:
    file "multiqc_data" into MULTIQC_DATA
    file "multiqc_report.html" into MULTIQC_REPORT

  when:
    params.output.multiqc == true

  script:
    """
    multiqc \
      --ignore ${workflow.launchDir}/${params.output.dir}/GEMs \
      --ignore ${workflow.launchDir}/${params.output.dir}/reports \
      ${workflow.launchDir}/${params.output.dir}
    """
}



/**
 * The createGEM process should run when all samples have
 * completed or if on a resume when the bootstrap signal is
 * received.
 */
CREATE_GEM_RUN = CREATE_GEM_READY_SIGNAL.mix(CREATE_GEM_BOOTSTRAP)

/**
 * Creates the GEM file from all the FPKM/TPM outputs
 */
process create_gem {
  label "python3"
  publishDir "${params.output.dir}/GEMs", mode: params.output.publish_mode

  input:
    val signal from CREATE_GEM_RUN.collect()

  output:
    file "*.GEM.*.txt" into GEM_FILES

  when:
    params.output.create_gem == true

  script:
  """
  create_gem.sh ${params.output.publish_fpkm} ${params.input.hisat2.enable} ${workflow.launchDir} ${params.output.dir} ${params.project.machine_name} ${params.output.publish_raw} ${params.output.publish_tpm}
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
    params.output.publish_sra == false

  script:
  """
  clean_work_files.sh ${files_list}
  """
}

/**
 * Merge the merged fastq files with the signals from hista2,
 * kallisto, salmon and fastqc_1 to clean up merged fastq files. This
 * is only needed for remote files that were downloaded
 * and then merged into a single sample in the fastq_merge
 * process.
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
    params.output.publish_downloaded_fastq == false

  script:
  """
  clean_work_files.sh ${files_list}
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
    params.output.publish_trimmed_fastq == false

  script:
  """
  clean_work_files.sh ${files_list}
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
    params.output.publish_sam == false

  script:
  """
  clean_work_files.sh ${files_list}
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
    params.output.publish_bam == false

  script:
  """
  clean_work_files.sh ${files_list}
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
    params.output.publish_gene_abundance == false

  script:
  """
  clean_work_dirs.sh ${directory}
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
    params.output.publish_gene_abundance == false

  script:
  """
  clean_work_files.sh ${files_list}
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
    params.output.publish_stringtie_gtf_and_ga == false

  script:
  """
  clean_work_files.sh ${files_list}
  """
}
