/*
========================================================================================
    VALIDATE INPUTS
========================================================================================
*/

def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

/**
 * Determine which quantification tool was selected.
 */
hisat2_enable = false
star_enable = false
kallisto_enable = false
salmon_enable = false

if (params.pipeline.equals('hisat2')) {
  hisat2_enable = true
}
else if (params.pipeline.equals('star')) {
  star_enable = true
}
else if (params.pipeline.equals('kallisto')) {
  kallisto_enable = true
}
else if (params.pipeline.equals('salmon')) {
  salmon_enable = true
}
else {
  error "Error: You must select a valid quantification tool using the '--pipeline' parameter. Currently valid options are 'salmon', 'kallisto', 'star', or 'hisat2'"
}

/**
 * Make sure that at least one output format is enabled.
 */
if ( (hisat2_enable && !params.hisat2_keep_counts && !params.hisat2_keep_fpkm && !params.hisat2_keep_tpm) ||
     (star_enable && !params.star_keep_counts && !params.star_keep_fpkm && !params.star_keep_tpm) ) {
  error "Error: at least one output format (raw, fpkm, tpm) must be enabled."
}

if (!hisat2_enable && !params.hisat2_keep_counts && !params.hisat2_keep_tpm) {
  error "Error: at least one output format (raw, tpm) must be enabled."
}

/**
 * Determine which output formats should be published.
 */
publish_fpkm = false
publish_tpm = false
publish_raw = false
publish_gem = false

if ((hisat2_enable && params.hisat2_keep_fpkm) ||
    (star_enable && params.star_keep_fpkm)) {
    publish_fpkm = true
}
if ((hisat2_enable && params.hisat2_keep_counts) ||
    (star_enable && params.star_keep_counts) ||
    (salmon_enable && params.salmon_keep_counts) ||
    (kallisto_enable && params.kallisto_keep_counts)) {
    publish_raw = true
}
if ((hisat2_enable && params.hisat2_keep_tpm) ||
    (star_enable && params.star_keep_tpm) ||
    (salmon_enable && params.salmon_keep_tpm) ||
    (kallisto_enable && params.kallisto_keep_tpm)) {
    publish_tpm = true
}
if ((hisat2_enable && params.hisat2_keep_gem) ||
    (star_enable && params.star_keep_gem) ||
    (salmon_enable && params.salmon_keep_gem) ||
    (kallisto_enable && params.kallisto_keep_gem)) {
    publish_gem = true
}

/**
 * Check that required reference files exist
 */
if (hisat2_enable && (!params.hisat2_gtf_file || file(params.hisat2_gtf_file).isEmpty())) {
    error "Error: GTF reference file for Hisat2 does not exist or is empty! Please Check that you have the proper references, that they are placed in the reference directory, and they are named properly.\
    \nGEMmaker is missing the following file: '${params.hisat2_gtf_file}' "
}

if (hisat2_enable && (!params.hisat2_index_dir || !file(params.hisat2_index_dir).isDirectory())) {
    error "Error: hisat2 Index Directory does not exist or is empty! Please Check that you have the proper references, that they are placed in the reference directory, and they are named properly.\
    \nGEMmaker is missing the following file: '${params.hisat2_index_dir}'"
}

if (hisat2_enable && !params.hisat2_base_name) {
    error "Error: hisat2 requires the hisat2 index file's base name! Please Check that you have assigned the base name parameter '--hisat2_base_name'"
}


if (star_enable && (!params.star_gtf_file || file(params.star_gtf_file).isEmpty())) {
    error "Error: GTF reference file for star does not exist or is empty! Please Check that you have the proper references, that they are placed in the reference directory, and they are named properly.\
    \nGEMmaker is missing the following file: '${params.star_gtf_file}' "
}

if (star_enable && (!params.star_index_dir || !file(params.star_index_dir).isDirectory())) {
    error "Error: star Index Directory does not exist or is empty! Please Check that you have the proper references, that they are placed in the reference directory, and they are named properly.\
    \nGEMmaker is missing the following file: '${params.star_index_dir}'"
}


if (kallisto_enable && (!params.kallisto_index_path || file(params.kallisto_index_path).isEmpty())) {
    error "Error: Kallisto Index File does not exist or is empty! Please Check that you have the proper references, that they are placed in the reference directory, and they are named properly.\
    \nGEMmaker is missing the following file: '${params.kallisto_index_path}'"
}


if (salmon_enable && (!params.salmon_index_path || !file(params.salmon_index_path).isDirectory())) {
    error "Error: Salmon Index Directory does not exist or is empty! Please Check that you have the proper references, that they are placed in the reference directory, and they are named properly.\
    \nGEMmaker is missing the following file: '${params.salmon_index_path}'"
}

/**
 * Check that other input files/directories exist
 */
if ((hisat2_enable || star_enable) && (!params.trimmomatic_clip_file || !file(params.trimmomatic_clip_file).exists())) {
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
 * Get list of samples to skip.
 */
skip_samples = params.skip_samples
    ? file(params.skip_samples).readLines().collect { line -> line.trim() }
    : []

// Validate input parameters
WorkflowGEMmaker.initialise(workflow, skip_samples, params, log)

// Check input path parameters to see if they exist
def checkPathParamList = [ params.input, params.sras]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Define the sentinel for "done" signals
DONE_SENTINEL = 1

/*
========================================================================================
    CONFIG FILES
========================================================================================
*/


/*
========================================================================================
    IMPORT LOCAL MODULES/SUBWORKFLOWS
========================================================================================
*/

// Don't overwrite global params.modules, create a copy instead and use that within the main script.
def modules = params.modules.clone()

//
// MODULE: Local to the pipeline
//
include { GET_SOFTWARE_VERSIONS as get_software_versions } from '../modules/local/get_software_versions' addParams( options: [publish_files : ['tsv':'']] )

// Module: clean_work_dirs
include { clean_work_dirs } from '../modules/local/clean_work_dirs' addParams()

// Module: clean_work_files
include { clean_work_files } from '../modules/local/clean_work_files' addParams()

// Module: create_gem
include { create_gem } from '../modules/local/create_gem' addParams(publish_fpkm: publish_fpkm, hisat2_enable: hisat2_enable, publish_tpm: publish_tpm, publish_raw: publish_raw)

// Module: download_runs
// Note: we have to pass in the work directory in pieces with a 'path:'
// prefix because Nextflow will recoginze the path, see that the
// diretory attributes have changed and won't use the cache. So, this is
// a bit of a hack to keep Nextflow from rerunning the retreive on a resume
// if it completed successfully.
include { download_runs } from '../modules/local/download_runs' addParams(workDirParent: "path:" + workflow.workDir.getParent(), workDirName: workflow.workDir.getName())

// Module: failed_run_report
include { failed_run_report } from '../modules/local/failed_run_report' addParams()

// Module: fastq_dump
publish_pattern_fastq_dump = params.keep_retrieved_fastq
    ? "{*.fastq}"
    : "{none}"
include { fastq_dump } from '../modules/local/fastq_dump' addParams(publish_pattern_fastq_dump: publish_pattern_fastq_dump, DONE_SENTINEL: DONE_SENTINEL)

// Module: fastqc
include { fastqc as fastqc_1 } from '../modules/local/fastqc' addParams(DONE_SENTINEL: DONE_SENTINEL)
include { fastqc as fastqc_2 } from '../modules/local/fastqc' addParams(DONE_SENTINEL: DONE_SENTINEL)

// Module: fastq_merge
include { fastq_merge } from '../modules/local/fastq_merge' addParams(DONE_SENTINEL: DONE_SENTINEL)

// Module: stringtie_fpkm_tpm
include { stringtie_fpkm_tpm } from '../modules/local/stringtie_fpkm_tpm' addParams(DONE_SENTINEL: DONE_SENTINEL,
    keep_fpkm: publish_fpkm, keep_tpm: publish_tpm, keep_counts: publish_raw)

// Module: hisat2
include { hisat2 } from '../modules/local/hisat2' addParams(DONE_SENTINEL: DONE_SENTINEL)

// Module: star
include { star } from '../modules/local/star' addParams(DONE_SENTINEL: DONE_SENTINEL)

// Module: kallisto
publish_pattern_kallisto_ga = params.kallisto_keep_data
    ? "{*.ga,*.log}"
    : "{*.log}"
include { kallisto } from '../modules/local/kallisto' addParams(publish_pattern_kallisto_ga: publish_pattern_kallisto_ga, DONE_SENTINEL: DONE_SENTINEL)

// Module: kallisto_tpm
include { kallisto_tpm } from '../modules/local/kallisto_tpm' addParams(DONE_SENTINEL: DONE_SENTINEL)

// Module: multiqc
include { multiqc } from '../modules/local/multiqc' addParams()

// Module: next_sample
include { next_sample } from '../modules/local/next_sample' addParams()

// Module: retrieve_sra_metadata
// Note: we have to pass in the work directory in pieces with a 'path:'
// prefix because Nextflow will recoginze the path, see that the
// diretory attributes have changed and won't use the cache. So, this is
// a bit of a hack to keep Nextflow from rerunning the retreive on a resume
// if it completed successfully.
include { retrieve_sra_metadata } from '../modules/local/retrieve_sra_metadata' addParams(workDirParent: "path:" + workflow.workDir.getParent(), workDirName: workflow.workDir.getName())

// Module: salmon
publish_pattern_salmon_ga = params.salmon_keep_data
    ? "{*.ga}"
    : "{*.ga/aux_info/meta_info.json,*.ga/libParams/flenDist.txt}"
include { salmon } from '../modules/local/salmon' addParams(publish_pattern_salmon_ga: publish_pattern_salmon_ga, DONE_SENTINEL: DONE_SENTINEL)

// Module: salmon_tpm
include { salmon_tpm } from '../modules/local/salmon_tpm' addParams(DONE_SENTINEL: DONE_SENTINEL)

// Module: samtools_index
// Hisat2
if (hisat2_enable) {
  publish_pattern_samtools_index = params.hisat2_keep_bam
      ? "{*.log,*.bam.bai}"
      : "{*.log}"
}
// STAR
else  {
  publish_pattern_samtools_index = params.star_keep_bam
      ? "{*.log,*.bam.bai}"
      : "{*.log}"
}
include { samtools_index } from '../modules/local/samtools_index' addParams(publish_pattern_samtools_index: publish_pattern_samtools_index)

// Module: samtools_sort
// Hisat2
if (hisat2_enable) {
  publish_pattern_samtools_sort = params.hisat2_keep_bam
      ? "{*.log,*.bam}"
      : "{*.log}"
}
// STAR
else  {
  publish_pattern_samtools_sort = params.star_keep_bam
      ? "{*.log,*.bam}"
      : "{*.log}"
}
include { samtools_sort } from '../modules/local/samtools_sort' addParams(publish_pattern_samtools_sort: publish_pattern_samtools_sort, DONE_SENTINEL: DONE_SENTINEL)

// Module: samtools_merge - used with STAR only
publish_pattern_samtools_merge = params.star_keep_bam
    ? "{*.log,*.bam}"
    : "{*.log}"
include { samtools_merge } from '../modules/local/samtools_merge' addParams(publish_pattern_samtools_merge: publish_pattern_samtools_merge, DONE_SENTINEL: DONE_SENTINEL)

// Module: stringtie
// Hisat2
if (hisat2_enable) {
  publish_pattern_stringtie_ga_gtf = params.hisat2_keep_data
      ? "{*.ga, *.gtf}"
      : "{none}"
}
// STAR
else  {
  publish_pattern_stringtie_ga_gtf = params.star_keep_data
      ? "{*.ga, *.gtf}"
      : "{none}"
}
include { stringtie } from '../modules/local/stringtie' addParams(publish_pattern_stringtie_ga_gtf: publish_pattern_stringtie_ga_gtf, DONE_SENTINEL: DONE_SENTINEL)

// Module: trimmomatic
publish_pattern_trimmomatic = params.trimmomatic_keep_trimmed_fastq
    ? "{*.trim.log,*_trim.fastq}"
    : "{*.trim.log}"
include { trimmomatic } from '../modules/local/trimmomatic' addParams(publish_pattern_trimmomatic: publish_pattern_trimmomatic, DONE_SENTINEL: DONE_SENTINEL)

/*
========================================================================================
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
========================================================================================
*/


/*
========================================================================================
    RUN MAIN WORKFLOW
========================================================================================
*/

workflow GEMmaker {
    CH_SOFTWARE_VERSIONS = Channel.empty()
    //
    // MODULE: Pipeline reporting
    //
    CH_SOFTWARE_VERSIONS
        .map { it -> if (it) [ it.baseName, it ] }
        .groupTuple()
        .map { it[1][0] }
        .flatten()
        .collect()
        .set { CH_SOFTWARE_VERSIONS }

    get_software_versions( CH_SOFTWARE_VERSIONS.map { it }.collect() )


    /**
     * Load local sample files.
     */
    LOCAL_SAMPLES = params.input
        ? Channel.fromFilePairs( params.input, size: -1 )
        : Channel.empty()

    LOCAL_SAMPLES_LIST = LOCAL_SAMPLES
        .map{ [it[0], it[1], "local"] }

    /**
     * Retrieve metadata for remote samples from NCBI SRA.
     */
    if ( params.sras ) {
        retrieve_sra_metadata(file(params.sras))

        /**
         * Parse remote samples from the SRR2SRX mapping.
         */
        REMOTE_SAMPLES_LIST = retrieve_sra_metadata.out.SRR2SRX
            .splitCsv()
            .groupTuple(by: 1)
            .map { [it[1], it[0].join(" "), "remote"] }
    }
    else {
        REMOTE_SAMPLES_LIST = Channel.empty()
    }

    /**
     * Prepare sample files for staging and processing.
     */
    LOCAL_SAMPLES_LIST
        .mix(REMOTE_SAMPLES_LIST)
        .filter { sample -> !skip_samples.contains(sample[0]) }
        .subscribe onNext: { sample ->

            // Stage samples that are not in the skip list.
            (sample_id, run_files_or_ids, sample_type) = sample

            // Don't create a sample file for samples that are already done.
            done_file = file("${workflow.workDir}/GEMmaker/done/${sample_id}.sample.csv")
            if (!done_file.exists()) {

                // If the sample file already exists don't add more to it.
                stage_file = file("${workflow.workDir}/GEMmaker/stage/${sample_id}.sample.csv")
                if (!stage_file.exists()) {
                    line = (sample_type == "local")
                        ? "${sample_id}\t${run_files_or_ids}\t${sample_type}"
                        : "${sample_id}\t${run_files_or_ids}\t${sample_type}"
                    stage_file << line
                }
            }
        }, onComplete: {

            // If we don't have any staged files this means we had no
            // samples or the workflow has been run already and it's done.
            num_staged = file("${workflow.workDir}/GEMmaker/stage/*.sample.csv").size()
            if (num_staged == 0) {
                done_file = file("${workflow.workDir}/GEMmaker/process/DONE")
                done_file << "DONE"
            }
            else {
                // Move the first batch of samples into processing.
                file("${workflow.workDir}/GEMmaker/stage/*.sample.csv")
                    .sort()
                    .take(params.max_cpus)
                    .each { sample ->
                        sample.moveTo("${workflow.workDir}/GEMmaker/process")
                    }
            }
        }

    /**
     * Watch the process directory for new samples. Each new
     * sample is sent to different downstream processes based
     * on whether it is a local or remote sample.
     *
     * Close the channel when the "done" file is received.
     */
    Channel.watchPath("${workflow.workDir}/GEMmaker/process")
        .until { it.name == "DONE" }
        .splitCsv(sep: '\t')
        .branch { sample ->
            local: sample[2] == "local"
            remote: sample[2] == "remote"
        }
        .set { NEXT_SAMPLE }

    /**
     * Merge local samples with local fastq files.
     */
    LOCAL_FASTQ_FILES = NEXT_SAMPLE.local
        .join(LOCAL_SAMPLES_LIST)
        .map { [it[0], it[3]] }

    /**
     * Download, extract, and merge remote samples.
     */
    if ( params.sras ) {
        REMOTE_SAMPLES = NEXT_SAMPLE.remote

        download_runs(REMOTE_SAMPLES)
        SRA_FILES = download_runs.out.SRA_FILES

        fastq_dump(SRA_FILES)
        DOWNLOADED_FASTQ_FILES = fastq_dump.out.FASTQ_FILES

        fastq_merge(DOWNLOADED_FASTQ_FILES)
        MERGED_FASTQ_FILES = fastq_merge.out.FASTQ_FILES

        FAILED_SAMPLES = Channel.empty()
            .mix(
                download_runs.out.FAILED_SAMPLES.map { it[0] },
                fastq_dump.out.FAILED_SAMPLES.map { it[0] })
    }
    else {
        REMOTE_SAMPLES = Channel.empty()
        MERGED_FASTQ_FILES = Channel.empty()
        FAILED_SAMPLES = Channel.empty()
    }

    /**
     * Combine local and remote samples.
     */
    FASTQ_FILES = LOCAL_FASTQ_FILES.mix(MERGED_FASTQ_FILES)

    /**
     * Perform fastqc analysis on all samples.
     */
    fastqc_1(FASTQ_FILES)

    /**
     * Process samples with hisat2 if enabled.
     */
    if ( hisat2_enable ) {
        // execute hisat2 pipeline
        HISAT2_INDEXES = Channel.fromPath( "${params.hisat2_index_dir}/*" ).collect()
        FASTA_ADAPTER = Channel.fromPath( params.trimmomatic_clip_file ).collect()
        GTF_FILE = Channel.fromPath( params.hisat2_gtf_file ).collect()

        trimmomatic(FASTQ_FILES, FASTA_ADAPTER)
        TRIMMED_FASTQ_FILES = trimmomatic.out.FASTQ_FILES
        MERGED_FASTQ_DONE = trimmomatic.out.DONE_SIGNAL

        fastqc_2(TRIMMED_FASTQ_FILES)

        hisat2(TRIMMED_FASTQ_FILES, HISAT2_INDEXES)
        SAM_FILES = hisat2.out.SAM_FILES

        samtools_sort(SAM_FILES)
        BAM_FILES = samtools_sort.out.BAM_FILES

        samtools_index(BAM_FILES)
        INDEXED_BAM_FILES = samtools_index.out.BAM_FILES

        stringtie(INDEXED_BAM_FILES, GTF_FILE)
        STRINGTIE_FILES = stringtie.out.GA_GTF_FILES

        stringtie_fpkm_tpm(STRINGTIE_FILES)
        RAW_FILES = stringtie_fpkm_tpm.out.RAW_FILES
        TPM_FILES = stringtie_fpkm_tpm.out.TPM_FILES
        FPKM_FILES = stringtie_fpkm_tpm.out.FPKM_FILES

        // prepare channels for downstream processes
        COMPLETED_SAMPLES = stringtie_fpkm_tpm.out.DONE_SIGNAL

        MULTIQC_FILES = Channel.empty()
            .mix(
                fastqc_1.out.REPORTS.flatten(),
                trimmomatic.out.LOGS.map { it[1] },
                fastqc_2.out.REPORTS.flatten(),
                hisat2.out.LOGS.map { it[1] },
                samtools_index.out.LOGS.map { it[1] })
    }

    /**
     * Process samples with star if enabled.
     */
    if ( star_enable ) {
        // execute star pipeline
        STAR_INDEXES = Channel.fromPath( params.star_index_dir ).collect()
        FASTA_ADAPTER = Channel.fromPath( params.trimmomatic_clip_file ).collect()
        GTF_FILE = Channel.fromPath( params.star_gtf_file ).collect()

        trimmomatic(FASTQ_FILES, FASTA_ADAPTER)
        TRIMMED_FASTQ_FILES = trimmomatic.out.FASTQ_FILES
        MERGED_FASTQ_DONE = trimmomatic.out.DONE_SIGNAL

        fastqc_2(TRIMMED_FASTQ_FILES)

        star(TRIMMED_FASTQ_FILES, STAR_INDEXES)
        SAM_FILES = star.out.SAM_FILES

        samtools_merge(SAM_FILES)
        MERGED_BAM_FILES = samtools_merge.out.BAM_FILES

        samtools_sort(MERGED_BAM_FILES)
        BAM_FILES = samtools_sort.out.BAM_FILES

        samtools_index(BAM_FILES)
        INDEXED_BAM_FILES = samtools_index.out.BAM_FILES

        stringtie(INDEXED_BAM_FILES, GTF_FILE)
        STRINGTIE_FILES = stringtie.out.GA_GTF_FILES

        stringtie_fpkm_tpm(STRINGTIE_FILES)
        RAW_FILES = stringtie_fpkm_tpm.out.RAW_FILES
        TPM_FILES = stringtie_fpkm_tpm.out.TPM_FILES
        FPKM_FILES = stringtie_fpkm_tpm.out.FPKM_FILES

        // prepare channels for downstream processes
        COMPLETED_SAMPLES = stringtie_fpkm_tpm.out.DONE_SIGNAL

        MULTIQC_FILES = Channel.empty()
            .mix(
                fastqc_1.out.REPORTS.flatten(),
                trimmomatic.out.LOGS.map { it[1] },
                fastqc_2.out.REPORTS.flatten(),
                star.out.LOGS.map { it[1] },
                samtools_index.out.LOGS.map { it[1] })
    }

    /**
     * Process samples with kallisto if enabled.
     */
    if ( kallisto_enable ) {
        // execute kallisto pipeline
        KALLISTO_INDEX = Channel.fromPath( params.kallisto_index_path ).collect()

        kallisto(FASTQ_FILES, KALLISTO_INDEX)
        GA_FILES = kallisto.out.GA_FILES
        MERGED_FASTQ_DONE = kallisto.out.DONE_SIGNAL

        kallisto_tpm(GA_FILES)
        RAW_FILES = kallisto_tpm.out.RAW_FILES
        TPM_FILES = kallisto_tpm.out.TPM_FILES
        FPKM_FILES = Channel.empty()

        // prepare channels for downstream processes
        COMPLETED_SAMPLES = kallisto_tpm.out.DONE_SIGNAL

        MULTIQC_FILES = Channel.empty()
            .mix(
                fastqc_1.out.REPORTS.flatten(),
                kallisto.out.LOGS.map { it[1] })
    }

    /**
     * Process samples with salmon if enabled.
     */
    if ( salmon_enable ) {
        // execute salmon pipeline
        SALMON_INDEXES = Channel.fromPath( params.salmon_index_path ).collect()

        salmon(FASTQ_FILES, SALMON_INDEXES)
        GA_FILES = salmon.out.GA_FILES
        MERGED_FASTQ_DONE = salmon.out.DONE_SIGNAL

        salmon_tpm(GA_FILES)
        RAW_FILES = salmon_tpm.out.RAW_FILES
        TPM_FILES = salmon_tpm.out.TPM_FILES
        FPKM_FILES = Channel.empty()

        // prepare channels for downstream processes
        COMPLETED_SAMPLES = salmon_tpm.out.DONE_SIGNAL

        MULTIQC_FILES = Channel.empty()
            .mix(
                fastqc_1.out.REPORTS.flatten(),
                salmon.out.LOGS.map { it[1] })
    }

    /**
     * Move a new sample into processing when a current
     * sample is completed or skipped.
     */
    SAMPLE_DONE_SIGNAL = Channel.empty()
        .mix(
            COMPLETED_SAMPLES.map { it[0] },
            FAILED_SAMPLES)

    next_sample(SAMPLE_DONE_SIGNAL)

    /**
     * Define a "bootstrap" signal for the post-processing steps
     * (multiqc, create_gem) to ensure that they will run even if
     * no new samples were processed via COMPLETED_SAMPLES, such as
     * when resuming a run in which all samples have already been
     * processed.
     */
    BOOTSTRAP_SIGNAL = Channel.fromList( [DONE_SENTINEL] )

    /**
     * The multiqc process should run when all samples have
     * completed or on a resume when the bootstrap signal is
     * received.
     */
    if ( params.publish_multiqc_report ) {
        MULTIQC_RUN = BOOTSTRAP_SIGNAL.mix(COMPLETED_SAMPLES)
        MULTIQC_CONFIG = Channel.fromPath( params.multiqc_config_file )
        MULTIQC_CUSTOM_LOGO = Channel.fromPath( params.multiqc_custom_logo )

        multiqc(
            MULTIQC_RUN.collect(),
            MULTIQC_CONFIG.collect(),
            MULTIQC_CUSTOM_LOGO.collect(),
            MULTIQC_FILES.collect())
    }

    /**
     * The create_gem process should run when all samples have
     * completed or if on a resume when the bootstrap signal is
     * received.
     */
    if ( publish_gem ) {
        CREATE_GEM_RUN = BOOTSTRAP_SIGNAL.mix(COMPLETED_SAMPLES)
        CREATE_GEM_FILES = Channel.empty()
            .mix(
                RAW_FILES,
                TPM_FILES,
                FPKM_FILES)

        create_gem(
            CREATE_GEM_RUN.collect(),
            CREATE_GEM_FILES.collect())
    }

    /**
     * Create a report of failed samples.
     */
    if ( params.sras ) {
        FAILED_RUNS = Channel.empty()
            .mix(
                retrieve_sra_metadata.out.FAILED_RUNS,
                download_runs.out.FAILED_RUNS.map { it[1] },
                fastq_dump.out.FAILED_RUNS.map { it[1] })

        FAILED_RUN_TEMPLATE = Channel.fromPath( params.failed_run_report_template )

        failed_run_report(
            FAILED_RUNS.collect(),
            FAILED_RUN_TEMPLATE.collect())
    }



    /**
     * CLEANING INTERMEDIATE FILES
     *
     * Intermediate files should be deleted once they are no
     * longer needed. However, if an intermediate file is an
     * input to a task, deleting the file will invalidate
     * Nextflow's cache and cause the task to be re-run, even
     * if the sample was already processed.
     *
     * To trick Nextflow, we truncate the file to zero size
     * while preserving the original size and read/write
     * timestamps in the file's metadata.
     *
     * We determine when to clean the files from a particular
     * channel by combining the channel with the "done" signals
     * of every process that uses it. For example, if a channel
     * is used by two downstream processes, we combine the source
     * channel with the two corresponding done signals, and then
     * use groupTuple(size: 3) to ensure that signals from all three
     * channels are received for a file before deleting it.
     *
     * The done signal should include the sample_id and the
     * DONE_SENTINEL to ensure that the total number of signals
     * for a file can be counted.
     */
    CLEAN_WORK_FILES = Channel.empty()
    CLEAN_WORK_DIRS = Channel.empty()

    /**
     * Clean sra files after they are used by fastq_dump.
     */
    if ( params.sras && !params.keep_sra ) {
        CLEAN_SRA_FILES = SRA_FILES
            .mix(fastq_dump.out.DONE_SIGNAL)
            .groupTuple(size: 2)
            .map { [it[0], it[1].flatten().findAll { v -> v != DONE_SENTINEL }] }

        CLEAN_WORK_FILES = CLEAN_WORK_FILES.mix(CLEAN_SRA_FILES)
    }

    /**
     * Clean downloaded fastq files after they are used by fastq_merge.
     */
    if ( params.sras && !params.keep_retrieved_fastq ) {
        CLEAN_DOWNLOADED_FASTQ_FILES = DOWNLOADED_FASTQ_FILES
            .mix(fastq_merge.out.DONE_SIGNAL)
            .groupTuple(size: 2)
            .map { [it[0], it[1].flatten().findAll { v -> v != DONE_SENTINEL }] }

        CLEAN_WORK_FILES = CLEAN_WORK_FILES.mix(CLEAN_DOWNLOADED_FASTQ_FILES)
    }

    /**
     * Clean merged fastq files after they are used by fastqc_1 and
     * the quantification tool.
     */
    if ( !params.keep_retrieved_fastq ) {
        CLEAN_MERGED_FASTQ_FILES = MERGED_FASTQ_FILES
            .mix(
                fastqc_1.out.DONE_SIGNAL,
                MERGED_FASTQ_DONE)
            .groupTuple(size: 3)
            .map { [it[0], it[1].flatten().findAll { v -> v != DONE_SENTINEL }] }

        CLEAN_WORK_FILES = CLEAN_WORK_FILES.mix(CLEAN_MERGED_FASTQ_FILES)
    }

    /**
     * Clean trimmed fastq files after they are used by fastqc_2 and hisat2
     */
    if ( hisat2_enable && !params.trimmomatic_keep_trimmed_fastq ) {
        CLEAN_TRIMMED_FASTQ_FILES = TRIMMED_FASTQ_FILES
            .mix(
                fastqc_2.out.DONE_SIGNAL,
                hisat2.out.DONE_SIGNAL)
            .groupTuple(size: 3)
            .map { [it[0], it[1].flatten().findAll { v -> v != DONE_SENTINEL }] }

        CLEAN_WORK_FILES = CLEAN_WORK_FILES.mix(CLEAN_TRIMMED_FASTQ_FILES)
    }

    /**
     * Clean trimmed fastq files after they are used by fastqc_2 and Star
     */
    if ( star_enable && !params.trimmomatic_keep_trimmed_fastq ) {
        CLEAN_TRIMMED_FASTQ_FILES = TRIMMED_FASTQ_FILES
            .mix(
                fastqc_2.out.DONE_SIGNAL,
                star.out.DONE_SIGNAL)
            .groupTuple(size: 3)
            .map { [it[0], it[1].flatten().findAll { v -> v != DONE_SENTINEL }] }

        CLEAN_WORK_FILES = CLEAN_WORK_FILES.mix(CLEAN_TRIMMED_FASTQ_FILES)
    }

    /**
     * Clean sam files after they are used by samtools_sort. Hisat and Star
     */
    if ( (hisat2_enable && !params.hisat2_keep_sam) ||
         (star_enable && !params.star_keep_sam) ) {
        CLEAN_SAM_FILES = SAM_FILES
            .mix(samtools_sort.out.DONE_SIGNAL)
            .groupTuple(size: 2)
            .map { [it[0], it[1].flatten().findAll { v -> v != DONE_SENTINEL }] }

        CLEAN_WORK_FILES = CLEAN_WORK_FILES.mix(CLEAN_SAM_FILES)
    }


    /**
     * Clean merged bam files after they are used by samtools_sort. Star
     */
    if ( star_enable && !params.star_keep_bam ) {
        CLEAN_MERGED_BAM_FILES = MERGED_BAM_FILES
            .mix(samtools_sort.out.DONE_SIGNAL)
            .groupTuple(size: 2)
            .map { [it[0], it[1].flatten().findAll { v -> v != DONE_SENTINEL }] }

        CLEAN_WORK_FILES = CLEAN_WORK_FILES.mix(CLEAN_MERGED_BAM_FILES)
    }


    /**
     * Clean bam files after they are used by stringtie.
     */
    if ( (hisat2_enable && !params.hisat2_keep_bam) ||
         (star_enable && !params.star_keep_bam) ) {
        CLEAN_BAM_FILES = BAM_FILES
            .mix(stringtie.out.DONE_SIGNAL)
            .groupTuple(size: 2)
            .map { [it[0], it[1].flatten().findAll { v -> v != DONE_SENTINEL }] }

        CLEAN_WORK_FILES = CLEAN_WORK_FILES.mix(CLEAN_BAM_FILES)
    }

    /**
     * Clean stringtie files after they are used by stringtie_fpkm_tpm.
     */
    if ( (hisat2_enable && !params.hisat2_keep_data) ||
         (star_enable && !params.star_keep_data) ) {
        CLEAN_STRINGTIE_FILES = STRINGTIE_FILES
            .mix(stringtie_fpkm_tpm.out.DONE_SIGNAL)
            .groupTuple(size: 2)
            .map { [it[0], it[1].flatten().findAll { v -> v != DONE_SENTINEL }] }

        CLEAN_WORK_FILES = CLEAN_WORK_FILES.mix(CLEAN_STRINGTIE_FILES)
    }

    /**
     * Clean kallisto files after they are used by kallisto_tpm.
     */
    if ( kallisto_enable && !params.kallisto_keep_data ) {
        CLEAN_KALLISTO_GA_FILES = GA_FILES
            .mix(kallisto_tpm.out.DONE_SIGNAL)
            .groupTuple(size: 2)
            .map { [it[0], it[1][0]] }

        CLEAN_WORK_DIRS = CLEAN_WORK_DIRS.mix(CLEAN_KALLISTO_GA_FILES)
    }

    /**
     * Clean salmon files after they are used by salmon_tpm.
     */
    if ( salmon_enable && !params.salmon_keep_data ) {
        CLEAN_SALMON_GA_FILES = GA_FILES
            .mix(salmon_tpm.out.DONE_SIGNAL)
            .groupTuple(size: 2)
            .map { [it[0], it[1].flatten().findAll { v -> v != DONE_SENTINEL }] }

        CLEAN_WORK_DIRS = CLEAN_WORK_DIRS.mix(CLEAN_SALMON_GA_FILES)
    }

    clean_work_files(CLEAN_WORK_FILES)
    clean_work_dirs(CLEAN_WORK_DIRS)
}
