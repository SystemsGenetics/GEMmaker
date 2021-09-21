#!/usr/bin/env nextflow

nextflow.enable.dsl=2
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

/**
 * Determine which quantification tool was selected.
 */
hisat2_enable = false
kallisto_enable = false
salmon_enable = false

if (params.pipeline.equals('hisat2')) {
  hisat2_enable = true
}
else if (params.pipeline.equals('kallisto')) {
  kallisto_enable = true
}
else if (params.pipeline.equals('salmon')) {
  salmon_enable = true
}
else {
  error "Error: You must select a valid quantification tool using the '--pipeline' parameter. Currently valid options are 'salmon', 'kallisto', or 'hisat2'"
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

/**
 * Determine which GEM formats should be published.
 */
publish_fpkm = false
publish_tpm = false
publish_raw = false
publish_gem = false

if (hisat2_enable && params.hisat2_keep_fpkm) {
    publish_fpkm = true
}
if ((hisat2_enable && params.hisat2_keep_counts) ||
    (salmon_enable && params.salmon_keep_counts) ||
    (kallisto_enable && params.kallisto_keep_counts)) {
    publish_raw = true
}
if ((hisat2_enable && params.hisat2_keep_tpm) ||
    (salmon_enable && params.salmon_keep_tpm) ||
    (kallisto_enable && params.kallisto_keep_tpm)) {
    publish_tpm = true
}
if ((hisat2_enable && params.hisat2_keep_gem) ||
    (salmon_enable && params.salmon_keep_gem) ||
    (kallisto_enable && params.kallisto_keep_gem)) {
    publish_gem = true
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
  Remote sample list path:    ${params.sras}
  Local sample glob:          ${params.input}
  Skip samples file:          ${params.skip_samples}


Reports
-------
  Report directory:           ${params.outdir}/reports


Quantification:
---------------
  Tool:                       ${params.pipeline}"""

if (hisat2_enable) {
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
else if (kallisto_enable) {
  println """\
  Kallisto Index File:        ${params.kallisto_index_path}
  """
}
else if (salmon_enable) {
  println """\
  Salmon Index Directory:     ${params.salmon_index_path}
  """
}

println """\

Published Results:
---------------
  Output Dir:                 ${params.outdir}/GEMs"""

if (publish_fpkm) {
    println """  FPKM counts:                Yes"""
}
if (publish_raw) {
    println """  Raw counts:                 Yes"""
}
if (publish_tpm) {
    println """  TPM counts:                 Yes"""
}
if (publish_gem) {
    println """  GEM file:                   Yes"""
}

println """\

"""



/**
 * Create the directories used for running batches
 */
file("${workflow.workDir}/GEMmaker").mkdir()
file("${workflow.workDir}/GEMmaker/stage").mkdir()
file("${workflow.workDir}/GEMmaker/process").mkdir()
file("${workflow.workDir}/GEMmaker/done").mkdir()

/**
 * Make sure that the user hasn't changed the quantification
 * tool on a resumed run.
 */
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
 * Check that required reference files exist
 */
if (hisat2_enable && file(params.hisat2_gtf_file).isEmpty()) {
    error "Error: GTF reference file for Hisat2 does not exist or is empty! Please Check that you have the proper references, that they are placed in the reference directory, and they are named properly.\
    \nGEMmaker is missing the following file: '${params.hisat2_gtf_file}' "
}

if (hisat2_enable && !file(params.hisat2_index_dir).isDirectory()) {
    error "Error: hisat2 Index Directory does not exist or is empty! Please Check that you have the proper references, that they are placed in the reference directory, and they are named properly.\
    \nGEMmaker is missing the following file: '${params.hisat2_index_dir}'"
}

if (kallisto_enable && file(params.kallisto_index_path).isEmpty()) {
    error "Error: Kallisto Index File does not exist or is empty! Please Check that you have the proper references, that they are placed in the reference directory, and they are named properly.\
    \nGEMmaker is missing the following file: '${params.kallisto_index_path}'"
}

if (salmon_enable && !file(params.salmon_index_path).isDirectory()) {
    error "Error: Salmon Index Directory does not exist or is empty! Please Check that you have the proper references, that they are placed in the reference directory, and they are named properly.\
    \nGEMmaker is missing the following file: '${params.salmon_index_path}'"
}

/**
 * Check that other input files/directories exist
 */
if (hisat2_enable && !file(params.trimmomatic_clip_file).exists()) {
    error "Error: The Trimmomatic clip file cannot be found at '${params.trimmomatic_clip_file}'."
}

if (params.sras && !file(params.sras).exists()) {
    error "Error: The NCBI download sample file does not exists at '${params.sras}'. This file must be provided. If you are not downloading samples from NCBI SRA the file must exist but can be left empty."
}

if (params.skip_samples && !file(params.skip_samples).exists()) {
    error "Error: The file containing samples to skip does not exists at '${params.skip_samples}'."
}

if (!file(params.failed_run_report_template).exists()) {
    error "Error: The failed run report template cannot be found at '${params.failed_run_report_template}'. This file comes with GEMmaker and is required."
}



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



workflow {
    get_software_versions()

    /**
     * Create value channels that can be reused
     */
    if (hisat2_enable) {
        HISAT2_INDEXES = Channel.fromPath("${params.hisat2_index_dir}/*").collect()
        FASTA_ADAPTER = Channel.fromPath("${params.trimmomatic_clip_file}").collect()
        GTF_FILE = Channel.fromPath("${params.hisat2_gtf_file}").collect()
    }
    else {
        HISAT2_INDEXES = Channel.empty()
        FASTA_ADAPTER = Channel.empty()
        GTF_FILE = Channel.empty()
    }

    if (kallisto_enable) {
        KALLISTO_INDEX = Channel.fromPath("${params.kallisto_index_path}").collect()
    }
    else {
        KALLISTO_INDEX = Channel.empty()
    }

    if (salmon_enable) {
        SALMON_INDEXES = Channel.fromPath("${params.salmon_index_path}").collect()
    }
    else {
        SALMON_INDEXES = Channel.empty()
    }

    FAILED_RUN_TEMPLATE = Channel.fromPath("${params.failed_run_report_template}").collect()
    MULTIQC_CONFIG = Channel.fromPath("${params.multiqc_config_file}").collect()
    MULTIQC_CUSTOM_LOGO = Channel.fromPath("${params.multiqc_custom_logo}").collect()

    if (params.skip_samples) {
        SKIP_SAMPLES_FILE = Channel.fromPath("${params.skip_samples}")
    }
    else {
        SKIP_SAMPLES_FILE = Channel.value('NA')
    }

    /**
     * Local Sample Input.
     * This checks the folder that the user has given
     */
    if (params.input) {
        LOCAL_SAMPLE_FILES_FOR_STAGING = Channel.fromFilePairs( "${params.input}", size: -1 )
        LOCAL_SAMPLE_FILES_FOR_JOIN = Channel.fromFilePairs( "${params.input}", size: -1 )
    }
    else {
        LOCAL_SAMPLE_FILES_FOR_STAGING = Channel.empty()
        LOCAL_SAMPLE_FILES_FOR_JOIN = Channel.empty()
    }

    /**
     * Remote fastq_run_id Input.
     */
    if (params.sras != "") {
        SRR_FILE = Channel.fromPath("${params.sras}")
    }
    else {
        SRR_FILE = Channel.empty()
    }

    retrieve_sra_metadata(SRR_FILE, SKIP_SAMPLES_FILE)
    REMOTE_SAMPLES_LIST = retrieve_sra_metadata.out.REMOTE_SAMPLES_LIST
    METADATA_FAILED_RUNS = retrieve_sra_metadata.out.METADATA_FAILED_RUNS

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

    // Clean up any files left over from a previous run by moving them
    // back to the stage directory.
    existing_files = file('work/GEMmaker/process/*')
    for (existing_file in existing_files) {
      existing_file.moveTo('work/GEMmaker/stage')
    }
    existing_files = null

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

    // If there are no staged files then the workflow will
    // end because it only proceeds when there are samples
    // in the processed directory.  However suppose the workflow
    // fails on multiqc and needs to be resumed.  The
    // following bootstraps the post-processsing portion of
    // the workflow
    if (staged_files.size() == 0) {
        MULTIQC_BOOTSTRAP = Channel.fromList( [1] )
        CREATE_GEM_BOOTSTRAP = Channel.fromList( [1] )
    }
    else {
        MULTIQC_BOOTSTRAP = Channel.empty()
        CREATE_GEM_BOOTSTRAP = Channel.empty()
    }
    staged_files = null

    write_stage_files(ALL_SAMPLES)
    SAMPLES_READY_SIGNAL = write_stage_files.out.SAMPLES_READY_SIGNAL

    // When all batch files are created we need to then
    // move the first file into the process directory.
    SAMPLES_READY_SIGNAL.collect().set { FIRST_SAMPLE_START_SIGNAL }

    // TODO: move exec code into workflow
    start_first_batch(FIRST_SAMPLE_START_SIGNAL)

    // Create the channel that will watch the process directory
    // for new files. When a new sample file is added
    // it will be read it and sent it through the workflow.
    NEXT_SAMPLE = Channel
       .watchPath("${workflow.workDir}/GEMmaker/process")

    read_sample_file(NEXT_SAMPLE)
    SAMPLE_FILE_CONTENTS = read_sample_file.out.SAMPLE_FILE_CONTENTS

    // Split our sample file contents into two different
    // channels, one for remote samples and another for local.
    SAMPLE_FILE_CONTENTS
      .splitCsv(quote: '"')
      .branch {
          local: it[2] =~ /local/
          remote: true
      }
      .set { SAMPLE_BRANCHES }

    LOCAL_SAMPLES = SAMPLE_BRANCHES.local
    REMOTE_SAMPLES = SAMPLE_BRANCHES.remote

    // Split our list of local samples into two pathways, onefor
    // FastQC analysis and the other for read counting.  We don't
    // do this for remote samples because they need downloading
    // first.
    LOCAL_SAMPLES = LOCAL_SAMPLES
      .map {[it[0], 'hi']}
      .mix(LOCAL_SAMPLE_FILES_FOR_JOIN)
      .groupTuple(size: 2)
      .map {[it[0], it[1][0]]}

    // Create the channels needed for signalling when
    // samples are completed.
    HISAT2_SAMPLE_COMPLETE_SIGNAL = Channel.empty()
    KALLISTO_SAMPLE_COMPLETE_SIGNAL = Channel.empty()
    SALMON_SAMPLE_COMPLETE_SIGNAL = Channel.empty()
    SKIP_DOWNLOAD_SAMPLE = Channel.empty()
    SKIP_DUMP_SAMPLE = Channel.empty()

    // Create the channel that will collate all the signals
    // and release a signal when the sample is complete
    SAMPLE_COMPLETE_SIGNAL = Channel.empty().mix(
        HISAT2_SAMPLE_COMPLETE_SIGNAL,
        KALLISTO_SAMPLE_COMPLETE_SIGNAL,
        SALMON_SAMPLE_COMPLETE_SIGNAL,
        SKIP_DUMP_SAMPLE.map { it[0] },
        SKIP_DOWNLOAD_SAMPLE.map { it[0] })

    // TODO: move exec code into workflow
    next_sample(SAMPLE_COMPLETE_SIGNAL)

    download_runs(REMOTE_SAMPLES)
    SRA_TO_EXTRACT = download_runs.out.SRA_TO_EXTRACT
    SRA_TO_CLEAN = download_runs.out.SRA_TO_CLEAN
    SKIP_DOWNLOAD_SAMPLE = download_runs.out.SKIP_DOWNLOAD_SAMPLE
    DOWNLOAD_FAILED_RUNS = download_runs.out.DOWNLOAD_FAILED_RUNS

    fastq_dump(SRA_TO_EXTRACT)
    DOWNLOADED_FASTQ_FOR_MERGING = fastq_dump.out.DOWNLOADED_FASTQ_FOR_MERGING
    DOWNLOADED_FASTQ_FOR_CLEANING = fastq_dump.out.DOWNLOADED_FASTQ_FOR_CLEANING
    SKIP_DUMP_SAMPLE = fastq_dump.out.SKIP_DUMP_SAMPLE
    CLEAN_SRA_SIGNAL = fastq_dump.out.CLEAN_SRA_SIGNAL
    FASTQ_DUMP_FAILED_RUNS = fastq_dump.out.FASTQ_DUMP_FAILED_RUNS

    fastq_merge(DOWNLOADED_FASTQ_FOR_MERGING)
    MERGED_SAMPLES_FOR_COUNTING = fastq_merge.out.MERGED_SAMPLES_FOR_COUNTING
    MERGED_SAMPLES_FOR_FASTQC_1 = fastq_merge.out.MERGED_SAMPLES_FOR_FASTQC_1
    MERGED_FASTQ_FOR_CLEANING = fastq_merge.out.MERGED_FASTQ_FOR_CLEANING
    CLEAN_FASTQ_SIGNAL = fastq_merge.out.CLEAN_FASTQ_SIGNAL

    /**
     * This is where we combine samples from both local and remote sources.
     */
    COMBINED_SAMPLES = LOCAL_SAMPLES.mix(MERGED_SAMPLES_FOR_FASTQC_1)

    fastqc_1(COMBINED_SAMPLES)
    FASTQC_1_OUTPUT = fastqc_1.out.FASTQC_1_OUTPUT
    CLEAN_MERGED_FASTQ_FASTQC_SIGNAL = fastqc_1.out.CLEAN_MERGED_FASTQ_FASTQC_SIGNAL

    /**
     * THIS IS WHERE THE SPLIT HAPPENS FOR hisat2 vs Kallisto vs Salmon
     *
     * Information about "choice" split operator (to be deleted before final
     * GEMmaker release)
     */
    COMBINED_SAMPLES
        .branch {
            hisat2:   hisat2_enable
            kallisto: kallisto_enable
            salmon:   salmon_enable
        }
        .set { SAMPLE_TOOL_BRANCHES }

    HISAT2_CHANNEL   = SAMPLE_TOOL_BRANCHES.hisat2
    KALLISTO_CHANNEL = SAMPLE_TOOL_BRANCHES.kallisto
    SALMON_CHANNEL   = SAMPLE_TOOL_BRANCHES.salmon

    kallisto(KALLISTO_CHANNEL, KALLISTO_INDEX)
    KALLISTO_GA = kallisto.out.KALLISTO_GA
    KALLISTO_GA_TO_CLEAN = kallisto.out.KALLISTO_GA_TO_CLEAN
    CLEAN_MERGED_FASTQ_KALLISTO_SIGNAL = kallisto.out.CLEAN_MERGED_FASTQ_KALLISTO_SIGNAL
    KALLISTO_LOG = kallisto.out.KALLISTO_LOG

    kallisto_tpm(KALLISTO_GA)
    KALLISTO_TPM = kallisto_tpm.out.KALLISTO_TPM
    KALLISTO_RAW = kallisto_tpm.out.KALLISTO_RAW
    CLEAN_KALLISTO_GA_SIGNAL = kallisto_tpm.out.CLEAN_KALLISTO_GA_SIGNAL
    KALLISTO_SAMPLE_COMPLETE_SIGNAL = kallisto_tpm.out.KALLISTO_SAMPLE_COMPLETE_SIGNAL

    salmon(SALMON_CHANNEL, SALMON_INDEXES)
    SALMON_GA = salmon.out.SALMON_GA
    SALMON_GA_LOG = salmon.out.SALMON_GA_LOG
    SALMON_GA_TO_CLEAN = salmon.out.SALMON_GA_TO_CLEAN
    CLEAN_MERGED_FASTQ_SALMON_SIGNAL = salmon.out.CLEAN_MERGED_FASTQ_SALMON_SIGNAL

    salmon_tpm(SALMON_GA)
    SALMON_TPM = salmon_tpm.out.SALMON_TPM
    SALMON_RAW = salmon_tpm.out.SALMON_RAW
    CLEAN_SALMON_GA_SIGNAL = salmon_tpm.out.CLEAN_SALMON_GA_SIGNAL
    SALMON_SAMPLE_COMPLETE_SIGNAL = salmon_tpm.out.SALMON_SAMPLE_COMPLETE_SIGNAL

    trimmomatic(HISAT2_CHANNEL, FASTA_ADAPTER)
    TRIMMED_SAMPLES_FOR_FASTQC = trimmomatic.out.TRIMMED_SAMPLES_FOR_FASTQC
    TRIMMED_SAMPLES_FOR_HISAT2 = trimmomatic.out.TRIMMED_SAMPLES_FOR_HISAT2
    TRIMMED_FASTQ_FOR_CLEANING = trimmomatic.out.TRIMMED_FASTQ_FOR_CLEANING
    TRIMMED_SAMPLE_LOG = trimmomatic.out.TRIMMED_SAMPLE_LOG

    fastqc_2(TRIMMED_SAMPLES_FOR_FASTQC)
    FASTQC_2_OUTPUT = fastqc_2.out.FASTQC_2_OUTPUT
    CLEAN_TRIMMED_FASTQ_FASTQC_SIGNAL = fastqc_2.out.CLEAN_TRIMMED_FASTQ_FASTQC_SIGNAL

    hisat2(TRIMMED_SAMPLES_FOR_HISAT2, HISAT2_INDEXES)
    HISAT2_SAM_FILE = hisat2.out.HISAT2_SAM_FILE
    HISAT2_SAM_LOG = hisat2.out.HISAT2_SAM_LOG
    SAM_FOR_CLEANING = hisat2.out.SAM_FOR_CLEANING
    CLEAN_TRIMMED_FASTQ_HISAT_SIGNAL = hisat2.out.CLEAN_TRIMMED_FASTQ_HISAT_SIGNAL
    CLEAN_MERGED_FASTQ_HISAT_SIGNAL = hisat2.out.CLEAN_MERGED_FASTQ_HISAT_SIGNAL

    samtools_sort(HISAT2_SAM_FILE)
    SORTED_FOR_INDEX = samtools_sort.out.SORTED_FOR_INDEX
    BAM_FOR_CLEANING = samtools_sort.out.BAM_FOR_CLEANING
    CLEAN_SAM_SIGNAL = samtools_sort.out.CLEAN_SAM_SIGNAL

    samtools_index(SORTED_FOR_INDEX)
    BAM_INDEXED_FOR_STRINGTIE = samtools_index.out.BAM_INDEXED_FOR_STRINGTIE
    BAI_INDEXED_FILE = samtools_index.out.BAI_INDEXED_FILE
    BAM_INDEXED_LOG = samtools_index.out.BAM_INDEXED_LOG

    stringtie(BAM_INDEXED_FOR_STRINGTIE, GTF_FILE)
    STRINGTIE_GTF_FOR_FPKM = stringtie.out.STRINGTIE_GTF_FOR_FPKM
    STRINGTIE_GTF_FOR_CLEANING = stringtie.out.STRINGTIE_GTF_FOR_CLEANING
    CLEAN_BAM_SIGNAL = stringtie.out.CLEAN_BAM_SIGNAL

    hisat2_fpkm_tpm(STRINGTIE_GTF_FOR_FPKM)
    HISAT2_FPKM = hisat2_fpkm_tpm.out.HISAT2_FPKM
    HISAT2_TPM = hisat2_fpkm_tpm.out.HISAT2_TPM
    HISAT2_RAW = hisat2_fpkm_tpm.out.HISAT2_RAW
    CLEAN_STRINGTIE_SIGNAL = hisat2_fpkm_tpm.out.CLEAN_STRINGTIE_SIGNAL
    HISAT2_SAMPLE_COMPLETE_SIGNAL = hisat2_fpkm_tpm.out.HISAT2_SAMPLE_COMPLETE_SIGNAL

    /**
     * The multiqc process should run when all samples have
     * completed or if on a resume when the bootstrap signal is
     * received.
     */
    MULTIQC_RUN = SAMPLE_COMPLETE_SIGNAL.mix(MULTIQC_BOOTSTRAP)

    FASTQC_1_OUTPUT.flatten()
      .concat(FASTQC_2_OUTPUT.flatten(), KALLISTO_LOG, SALMON_GA_LOG, HISAT2_SAM_LOG, BAM_INDEXED_LOG, TRIMMED_SAMPLE_LOG).set {MULTIQC_FILES}

    if ( params.publish_multiqc_report == true ) {
        multiqc(MULTIQC_RUN.collect(), MULTIQC_CONFIG, MULTIQC_CUSTOM_LOGO, MULTIQC_FILES.collect())
        MULTIQC_DATA = multiqc.out.MULTIQC_DATA
        MULTIQC_REPORT = multiqc.out.MULTIQC_REPORT
    }

    /**
     * The create_gem process should run when all samples have
     * completed or if on a resume when the bootstrap signal is
     * received.
     */
    CREATE_GEM_RUN = SAMPLE_COMPLETE_SIGNAL.mix(CREATE_GEM_BOOTSTRAP)
    SALMON_RAW
      .concat(SALMON_TPM, KALLISTO_RAW, KALLISTO_TPM, HISAT2_RAW, HISAT2_TPM, HISAT2_FPKM).set {QUANT_FILES}

    if ( publish_gem == true ) {
        create_gem(CREATE_GEM_RUN.collect(), QUANT_FILES.collect())
        GEM_FILES = create_gem.out.GEM_FILES
    }

    failed_run_report(METADATA_FAILED_RUNS, DOWNLOAD_FAILED_RUNS.collect(), FASTQ_DUMP_FAILED_RUNS.collect(), FAILED_RUN_TEMPLATE)
    FAILED_RUN_REPORT = failed_run_report.out.FAILED_RUN_REPORT

    SRA_TO_CLEAN
      .mix(CLEAN_SRA_SIGNAL)
      .groupTuple(size: 2)
      .set { CLEAN_SRA_READY }

    if ( params.keep_sra == false ) {
        clean_sra(CLEAN_SRA_READY)
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

    FCLEAN = DOWNLOADED_FASTQ_FOR_CLEANING.mix(CLEAN_FASTQ_SIGNAL)
    FCLEAN.groupTuple(size: 2).set { FASTQ_CLEANUP_READY }

    if ( params.keep_retrieved_fastq == false ) {
        clean_fastq(FASTQ_CLEANUP_READY)
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

    if ( params.keep_retrieved_fastq == false ) {
        clean_merged_fastq(MERGED_FASTQ_CLEANUP_READY)
    }

    /**
     * Merge the Trimmomatic samples with Hisat's signal and FastQC signal. Once
     * both tools send the signal that they are done with the trimmed file it can
     * be removed.
     */
    TRHIMIX = TRIMMED_FASTQ_FOR_CLEANING.mix(CLEAN_TRIMMED_FASTQ_HISAT_SIGNAL, CLEAN_TRIMMED_FASTQ_FASTQC_SIGNAL)
    TRHIMIX.groupTuple(size: 3).set { TRIMMED_FASTQ_CLEANUP_READY }

    if ( params.trimmomatic_keep_trimmed_fastq == false ) {
        clean_trimmed_fastq(TRIMMED_FASTQ_CLEANUP_READY)
    }

    /**
     * Merge the HISAT sam file with samtools_sort signal that it is
     * done so that we can remove these files.
     */
    HISSMIX = SAM_FOR_CLEANING.mix(CLEAN_SAM_SIGNAL)
    HISSMIX.groupTuple(size: 2).set { SAM_CLEANUP_READY }

    if ( params.hisat2_keep_sam == false ) {
        clean_sam(SAM_CLEANUP_READY)
    }

    /**
     * Merge the samtools_sort bam file with stringtie signal that it is
     * done so that we can remove these files.
     */
    SSSTMIX = BAM_FOR_CLEANING.mix(CLEAN_BAM_SIGNAL)
    SSSTMIX.groupTuple(size: 2).set { BAM_CLEANUP_READY }

    if ( params.hisat2_keep_bam == false ) {
        clean_bam(BAM_CLEANUP_READY)
    }

    /**
     * Merge the Kallisto .ga file with the clean kallisto ga signal so that we can
     * remove the .ga file after it has been used
     */
    KGAMIX = KALLISTO_GA_TO_CLEAN.mix(CLEAN_KALLISTO_GA_SIGNAL)
    KGAMIX.groupTuple(size: 2).set { KALLISTO_GA_CLEANUP_READY }

    if ( params.kallisto_keep_data == false ) {
        clean_kallisto_ga(KALLISTO_GA_CLEANUP_READY)
    }

    /**
     * Merge the Salmon .ga file with the clean salmon ga signal so that we can
     * remove the .ga file after it has been used
     */
    SGAMIX = SALMON_GA_TO_CLEAN.mix(CLEAN_SALMON_GA_SIGNAL)
    SGAMIX.groupTuple(size: 2).set { SALMON_GA_CLEANUP_READY }

    if ( params.salmon_keep_data == false ) {
        clean_salmon_ga(SALMON_GA_CLEANUP_READY)
    }

    /**
     * Merge the Salmon .ga file with the clean salmon ga signal so that we can
     * remove the .ga file after it has been used
     */
    SMIX = STRINGTIE_GTF_FOR_CLEANING.mix(CLEAN_STRINGTIE_SIGNAL)
    SMIX.groupTuple(size: 2).set { STRINGTIE_CLEANUP_READY }

    if ( params.hisat2_keep_data == false ) {
        clean_stringtie_ga(STRINGTIE_CLEANUP_READY)
    }

}



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
    path 'software_versions_mqc.yaml', emit: ch_software_versions_yaml
    path 'software_versions.csv'

    script:
    """
    echo $workflow.manifest.version > v_pipeline.txt
    echo $workflow.nextflow.version > v_nextflow.txt
    fastqc --version > v_fastqc.txt
    multiqc --version > v_multiqc.txt
    kallisto version > v_kallisto.txt
    hisat2 --version > v_hisat2.txt
    salmon --version > v_salmon.txt
    python --version > v_python.txt
    samtools version > v_samtools.txt
    fastq-dump --version > v_fastq_dump.txt
    stringtie --version > v_stringtie.txt
    trimmomatic -version > v_trimmomatic.txt
    pip freeze | grep pandas > v_pandas.txt
    scrape_software_versions.py &> software_versions_mqc.yaml
    """
}



/**
 * Retrieves metadata for all of the remote samples
 * and maps SRA runs to SRA experiments.
 */
process retrieve_sra_metadata {
  publishDir params.outdir, mode: params.publish_dir_mode, pattern: "failed_runs.metadata.txt"
  label "retrieve_sra_metadata"

  input:
    path(srr_file)
    path(skip_samples)

  output:
    stdout emit: REMOTE_SAMPLES_LIST
    path("failed_runs.metadata.txt"), emit: METADATA_FAILED_RUNS

  script:
  """
  >&2 echo "#TRACE n_remote_run_ids=`cat ${srr_file} | wc -l`"

  retrieve_sra_metadata.py \
      --run_id_file ${srr_file} \
      --meta_dir ${workflow.workDir}/GEMmaker \
      ${skip_samples != 'NA' ? "--skip_file ${skip_samples}" : ""}
  """
}



/**
 * Writes the batch files and stores them in the
 * stage directory.
 */
process write_stage_files {
  tag { sample_id }
  label "local"

  input:
    tuple val(sample_id), val(run_files_or_ids), val(sample_type)

  output:
    val(1), emit: SAMPLES_READY_SIGNAL

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



/**
 * Moves the first set of sample files into the process directory.
 */
process start_first_batch {
  label "local"
  cache false

  input:
    val(signal)

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
      // NEXT_SAMPLE.close()
      // HISAT2_SAMPLE_COMPLETE_SIGNAL.close()
      // KALLISTO_SAMPLE_COMPLETE_SIGNAL.close()
      // SALMON_SAMPLE_COMPLETE_SIGNAL.close()
      // SAMPLE_COMPLETE_SIGNAL.close()
      // MULTIQC_BOOTSTRAP.close()
      // CREATE_GEM_BOOTSTRAP.close()
      // SKIP_DUMP_SAMPLE.close()
      // SKIP_DOWNLOAD_SAMPLE.close()
      println "There are no staged samples.  Moving on to post-processing"
   }
}



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
    path(sample_file)

  output:
    stdout emit: SAMPLE_FILE_CONTENTS

  script:
  """
  cat ${sample_file}
  """
}



/**
 * Handles the end of a sample by moving a new sample
 * file into the process directory which triggers
 * the NEXT_SAMPLE.watchPath channel.
 */
process next_sample {
  tag { sample_id }
  label "local"
  cache false
  maxForks 1

  input:
    val(sample_id)

  exec:
    // Move the completed file into the done folder.
    sample_file = file("${workflow.workDir}/GEMmaker/process/" + sample_id + '.sample.csv')
    sample_file.moveTo("${workflow.workDir}/GEMmaker/done")

    // Move the next sample file into the processing directory
    // which will trigger the start of the next sample.
    staged_files = file("${workflow.workDir}/GEMmaker/stage/*")
    if (staged_files.size() > 0) {
        staged_files.first().moveTo("${workflow.workDir}/GEMmaker/process")
    }
    // If there are no more staged files, then check if there are no more
    // files in processing. If so, then close out all the channels.
    else {
        processing_files = file("${workflow.workDir}/GEMmaker/process/*.sample.csv")
        if (processing_files.size() == 0) {
            // NEXT_SAMPLE.close()
            // SAMPLE_COMPLETE_SIGNAL.close()
            // HISAT2_SAMPLE_COMPLETE_SIGNAL.close()
            // KALLISTO_SAMPLE_COMPLETE_SIGNAL.close()
            // SALMON_SAMPLE_COMPLETE_SIGNAL.close()
            // SAMPLE_COMPLETE_SIGNAL.close()
            // MULTIQC_BOOTSTRAP.close()
            // CREATE_GEM_BOOTSTRAP.close()
            // SKIP_DUMP_SAMPLE.close()
            // SKIP_DOWNLOAD_SAMPLE.close()
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
    tuple val(sample_id), val(run_ids), val(type)

  output:
    tuple val(sample_id), path("*.sra"), optional: true, emit: SRA_TO_EXTRACT
    tuple val(sample_id), path("*.sra"), optional: true, emit: SRA_TO_CLEAN
    tuple val(sample_id), path("sample_failed"), optional: true, emit: SKIP_DOWNLOAD_SAMPLE
    tuple val(sample_id), path('*.failed_runs.download.txt'), emit: DOWNLOAD_FAILED_RUNS

  script:
  """
  echo "#TRACE sample_id=${sample_id}"
  echo "#TRACE n_remote_run_ids=${run_ids.tokenize(' ').size()}"
  echo "#TRACE n_spots=`retrieve_sra_spots.py ${workflow.workDir}/GEMmaker ${sample_id}`"

  retrieve_sra.py --sample ${sample_id} --run_ids ${run_ids} --akey \${ASPERA_KEY}
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
    tuple val(sample_id), path(sra_files)

  output:
    tuple val(sample_id), path("*.fastq"), optional: true, emit: DOWNLOADED_FASTQ_FOR_MERGING
    tuple val(sample_id), path("*.fastq"), optional: true, emit: DOWNLOADED_FASTQ_FOR_CLEANING
    tuple val(sample_id), path("sample_failed"), optional: true, emit: SKIP_DUMP_SAMPLE
    tuple val(sample_id), val(1), emit: CLEAN_SRA_SIGNAL
    tuple val(sample_id), path('*.failed_runs.fastq-dump.txt'), emit: FASTQ_DUMP_FAILED_RUNS

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
    tuple val(sample_id), path(fastq_files)

  output:
    tuple val(sample_id), path("${sample_id}_?.fastq"), emit: MERGED_SAMPLES_FOR_COUNTING
    tuple val(sample_id), path("${sample_id}_?.fastq"), emit: MERGED_SAMPLES_FOR_FASTQC_1
    tuple val(sample_id), path("${sample_id}_?.fastq"), emit: MERGED_FASTQ_FOR_CLEANING
    tuple val(sample_id), val(1), emit: CLEAN_FASTQ_SIGNAL

  script:
  """
  echo "#TRACE sample_id=${sample_id}"
  echo "#TRACE fastq_lines=`cat *.fastq | wc -l`"

  fastq_merge.sh ${sample_id}
  """
}



/**
 * Performs fastqc on raw fastq files
 */
process fastqc_1 {
  publishDir "${params.outdir}/Samples/${sample_id}", mode: params.publish_dir_mode, pattern: "*_fastqc.*"
  tag { sample_id }
  label "fastqc"

  input:
    tuple val(sample_id), path(fastq_files)

  output:
    tuple path("*_fastqc.html"), path("*_fastqc.zip"), emit: FASTQC_1_OUTPUT
    tuple val(sample_id), val(1), emit: CLEAN_MERGED_FASTQ_FASTQC_SIGNAL

  script:
  """
  echo "#TRACE sample_id=${sample_id}"
  echo "#TRACE fastq_lines=`cat *.fastq | wc -l`"

  fastqc ${fastq_files}
  """
}



/**
 * Performs KALLISTO alignemnt of fastq files
 */
process kallisto {
  publishDir "${params.outdir}/Samples/${sample_id}", mode: params.publish_dir_mode, pattern: publish_pattern_Kallisto_GA
  tag { sample_id }
  label "process_medium"
  label "kallisto"

  input:
    tuple val(sample_id), path(fastq_files)
    path(kallisto_index)

  output:
    tuple val(sample_id), path("*.ga"), emit: KALLISTO_GA
    tuple val(sample_id), path("*.ga"), emit: KALLISTO_GA_TO_CLEAN
    tuple val(sample_id), val(1), emit: CLEAN_MERGED_FASTQ_KALLISTO_SIGNAL
    path("*.kallisto.log"), emit: KALLISTO_LOG

  script:
  """
  echo "#TRACE sample_id=${sample_id}"
  echo "#TRACE fastq_lines=`cat *.fastq | wc -l`"
  echo "#TRACE index_bytes=`stat -Lc '%s' ${kallisto_index}`"

  kallisto.sh \
    ${sample_id} \
    ${kallisto_index} \
    ${params.kallisto_bootstrap_samples} \
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
    tuple val(sample_id), path(ga_File)

  output:
    path("*.tpm"), optional: true, emit: KALLISTO_TPM
    path("*.raw"), optional: true, emit: KALLISTO_RAW
    tuple val(sample_id), val(1), emit: CLEAN_KALLISTO_GA_SIGNAL
    val(sample_id), emit: KALLISTO_SAMPLE_COMPLETE_SIGNAL

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
  label "process_medium"
  label "salmon"

  input:
    tuple val(sample_id), path(fastq_files)
    path(salmon_index)

  output:
    tuple val(sample_id), path("*.ga"), emit: SALMON_GA
    tuple val(sample_id), path("${sample_id}.Salmon.ga", type: 'dir'), emit: SALMON_GA_LOG
    tuple val(sample_id), path("*.ga/quant.sf"), emit: SALMON_GA_TO_CLEAN
    tuple val(sample_id), val(1), emit: CLEAN_MERGED_FASTQ_SALMON_SIGNAL

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
    tuple val(sample_id), path(ga_file)

  output:
    path("*.Salmon.tpm"), optional: true, emit: SALMON_TPM
    path("*.Salmon.raw"), optional: true, emit: SALMON_RAW
    tuple val(sample_id), val(1), emit: CLEAN_SALMON_GA_SIGNAL
    val(sample_id), emit: SALMON_SAMPLE_COMPLETE_SIGNAL

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
  label "process_medium"
  label "trimmomatic"

  input:
    tuple val(sample_id), path(fastq_files)
    path(fasta_adapter)

  output:
    tuple val(sample_id), path("*_trim.fastq"), emit: TRIMMED_SAMPLES_FOR_FASTQC
    tuple val(sample_id), path("*_trim.fastq"), emit: TRIMMED_SAMPLES_FOR_HISAT2
    tuple val(sample_id), path("*_trim.fastq"), emit: TRIMMED_FASTQ_FOR_CLEANING
    tuple val(sample_id), path("*.trim.log"), emit: TRIMMED_SAMPLE_LOG

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
    tuple val(sample_id), path(fastq_files)

  output:
    tuple path("*_fastqc.html"), path("*_fastqc.zip"), emit: FASTQC_2_OUTPUT
    tuple val(sample_id), val(1), emit: CLEAN_TRIMMED_FASTQ_FASTQC_SIGNAL

  script:
  """
  echo "#TRACE sample_id=${sample_id}"
  echo "#TRACE trimmed_fastq_lines=`cat *.fastq | wc -l`"

  fastqc ${fastq_files}
  """
}



/**
 * Performs hisat2 alignment of fastq files to a genome reference
 */
process hisat2 {
  publishDir "${params.outdir}/Samples/${sample_id}", mode: params.publish_dir_mode, pattern: "*.log"
  tag { sample_id }
  label "process_medium"
  label "hisat2"

  input:
    tuple val(sample_id), path(fastq_files)
    path(indexes)

  output:
    tuple val(sample_id), path("*.sam"), emit: HISAT2_SAM_FILE
    tuple val(sample_id), path("*.sam.log"), emit: HISAT2_SAM_LOG
    tuple val(sample_id), path("*.sam"), emit: SAM_FOR_CLEANING
    tuple val(sample_id), val(1), emit: CLEAN_TRIMMED_FASTQ_HISAT_SIGNAL
    tuple val(sample_id), val(1), emit: CLEAN_MERGED_FASTQ_HISAT_SIGNAL

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
 */
process samtools_sort {
  publishDir "${params.outdir}/Samples/${sample_id}", mode: params.publish_dir_mode, pattern: publish_pattern_samtools_sort
  tag { sample_id }
  label "samtools"

  input:
    tuple val(sample_id), path(sam_file)

  output:
    tuple val(sample_id), path("*.bam"), emit: SORTED_FOR_INDEX
    tuple val(sample_id), path("*.bam"), emit: BAM_FOR_CLEANING
    tuple val(sample_id), val(1), emit: CLEAN_SAM_SIGNAL

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
 */
process samtools_index {
  publishDir "${params.outdir}/Samples/${sample_id}", mode: params.publish_dir_mode, pattern: publish_pattern_samtools_index
  tag { sample_id }
  label "samtools"

  input:
    tuple val(sample_id), path(bam_file)

  output:
    tuple val(sample_id), path(bam_file), emit: BAM_INDEXED_FOR_STRINGTIE
    tuple val(sample_id), path("*.bam.bai"), emit: BAI_INDEXED_FILE
    tuple val(sample_id), path("*.bam.log"), emit: BAM_INDEXED_LOG

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
 */
process stringtie {
  publishDir "${params.outdir}/Samples/${sample_id}", mode: params.publish_dir_mode, pattern: publish_pattern_stringtie_gtf_and_ga
  tag { sample_id }
  label "process_medium"
  label "stringtie"

  input:
    tuple val(sample_id), path(bam_file)
    path(gtf_file)

  output:
    tuple val(sample_id), path("*.Hisat2.ga"), path("*.Hisat2.gtf"), emit: STRINGTIE_GTF_FOR_FPKM
    tuple val(sample_id), path("*.Hisat2.*"), emit: STRINGTIE_GTF_FOR_CLEANING
    tuple val(sample_id), val(1), emit: CLEAN_BAM_SIGNAL

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
    tuple val(sample_id), path(ga_file), path(gtf_file)

  output:
    path("*.Hisat2.fpkm"), optional: true, emit: HISAT2_FPKM
    path("*.Hisat2.tpm"), optional: true, emit: HISAT2_TPM
    path("*.Hisat2.raw"), optional: true, emit: HISAT2_RAW
    tuple val(sample_id), val(1), emit: CLEAN_STRINGTIE_SIGNAL
    val(sample_id), emit: HISAT2_SAMPLE_COMPLETE_SIGNAL

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
 * Process to generate the multiqc report once everything is completed
 */
process multiqc {
  label "multiqc"
  cache false
  publishDir "${params.outdir}/reports", mode: params.publish_dir_mode

  input:
    val(signal)
    path(multiqc_config)
    path(gemmaker_logo)
    path(input_files)

  output:
    path("multiqc_data"), emit: MULTIQC_DATA
    path("multiqc_report.html"), emit: MULTIQC_REPORT

  script:
  """
  multiqc --config ${multiqc_config} ./
  """
}



/**
 * Creates the GEM file from all the FPKM/TPM outputs
 */
process create_gem {
  label "create_gem"
  publishDir "${params.outdir}/GEMs", mode: params.publish_dir_mode

  input:
    val(signal)
    path(input_files)

  output:
    path("*.GEM.*.txt"), emit: GEM_FILES

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
    path(metadata_failed_runs)
    path(download_failed_runs)
    path(fastq_dump_failed_runs)
    path(failed_run_template)

  output:
    path("failed_SRA_run_report.html"), emit: FAILED_RUN_REPORT

  script:
  """
  failed_runs_report.py --template ${failed_run_template}
  """
}



/**
 * Cleans downloaded SRA files
 */
process clean_sra {
  tag { sample_id }

  input:
    tuple val(sample_id), val(files_list)

  script:
  """
  clean_work_files.sh "${files_list[0]}"
  """
}



/**
 * Cleans downloaded fastq files
 */
process clean_fastq {
  tag { sample_id }

  input:
    tuple val(sample_id), val(files_list)

  script:
  """
  clean_work_files.sh "${files_list[0].join(" ")}"
  """
}



/**
 * Cleans merged fastq files
 */
process clean_merged_fastq {
  tag { sample_id }

  input:
    tuple val(sample_id), val(files_list)

  script:
  """
  clean_work_files.sh "${files_list[0].join(" ")}"
  """
}



/**
 * Cleans trimmed fastq files
 */
process clean_trimmed_fastq {
  tag { sample_id }

  input:
    tuple val(sample_id), val(files_list)

  script:
  """
  clean_work_files.sh "${files_list[0].join(" ")}"
  """
}



/**
 * Clean up SAM files
 */
process clean_sam {
  tag { sample_id }

  input:
    tuple val(sample_id), val(files_list)

  script:
  """
  clean_work_files.sh "${files_list[0]}"
  """
}



/**
 * Clean up BAM files
 */
process clean_bam {
  tag { sample_id }

  input:
    tuple val(sample_id), val(files_list)

  script:
  """
  clean_work_files.sh "${files_list[0]}"
  """
}



/**
 * Clean up Kallisto GA files
 */
process clean_kallisto_ga {
  tag { sample_id }

  input:
    tuple val(sample_id), val(directory)

  script:
  """
  clean_work_dirs.sh "${directory[0]}"
  """
}



/**
 * Clean up Salmon GA files
 */
process clean_salmon_ga {
  tag { sample_id }

  input:
    tuple val(sample_id), val(files_list)

  script:
  """
  clean_work_files.sh "${files_list[0]}"
  """
}



/**
 * Clean up Salmon GA files
 */
process clean_stringtie_ga {
  tag { sample_id }

  input:
    tuple val(sample_id), val(files_list)

  script:
  """
  clean_work_files.sh "${files_list[0].join(" ")}"
  """
}
