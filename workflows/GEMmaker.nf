/*
========================================================================================
    VALIDATE INPUTS
========================================================================================
*/

def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

// Validate input parameters
WorkflowGEMmaker.initialise(workflow, params, log)

// Check input path parameters to see if they exist
def checkPathParamList = [ params.input ]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

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
if (hisat2_enable && !params.hisat2_keep_counts && !params.hisat2_keep_fpkm && !params.hisat2_keep_tpm) {
  error "Error: at least one output format (raw, fpkm, tpm) must be enabled for hisat2"
}

if (!hisat2_enable && !params.hisat2_keep_counts && !params.hisat2_keep_tpm) {
  error "Error: at least one output format (raw, tpm) must be enabled for kallisto / salmon"
}

/**
 * Determine which output formats should be published.
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
include { GET_SOFTWARE_VERSIONS } from '../modules/local/get_software_versions' addParams( options: [publish_files : ['tsv':'']] )

/*
========================================================================================
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
========================================================================================
*/


/*
========================================================================================
    PRE-WORKFLOW SETUP
========================================================================================
*/

/**
 * Define publish pattern for downloaded FASTQ files
 */
publish_pattern_fastq_dump = params.keep_retrieved_fastq
  ? "{*.fastq}"
  : "{none}"

/**
 * Define publish pattern for trimmed FASTQ files
 */
publish_pattern_trimmomatic = params.trimmomatic_keep_trimmed_fastq
  ? "{*.trim.log,*_trim.fastq}"
  : "{*.trim.log}"

/**
 * Define publish pattern for BAM files
 */
publish_pattern_samtools_sort = params.hisat2_keep_bam
  ? "{*.log,*.bam}"
  : "{*.log}"

publish_pattern_samtools_index = params.hisat2_keep_bam
  ? "{*.log,*.bam.bai}"
  : "{*.log}"

/**
 * Define publish pattern for stringtie GA and GTF files
 */
publish_pattern_stringtie_ga_gtf = params.hisat2_keep_data
  ? "{*.ga, *.gtf}"
  : "{none}"

/**
 * Define publish pattern for Kallisto GA files
 */
publish_pattern_kallisto_ga = params.kallisto_keep_data
  ? "{*.ga,*.log}"
  : "{*.log}"

/**
 * Define publish pattern for Salmon GA files
 * Publishes only log file used by multiqc if false
 */
publish_pattern_salmon_ga = params.salmon_keep_data
  ? "{*.ga}"
  : "{*.ga/aux_info/meta_info.json,*.ga/libParams/flenDist.txt}"


/**
 * Define the sentinel for "done" signals
 */
DONE_SENTINEL = 1


/*
========================================================================================
    RUN MAIN WORKFLOW
========================================================================================
*/

workflow GEMmaker {
    ch_software_versions = Channel.empty()


    //
    // MODULE: Pipeline reporting
    //
    ch_software_versions
        .map { it -> if (it) [ it.baseName, it ] }
        .groupTuple()
        .map { it[1][0] }
        .flatten()
        .collect()
        .set { ch_software_versions }

    GET_SOFTWARE_VERSIONS (
        ch_software_versions.map { it }.collect()
    )

    /**
     * Delete "done" file from previous run.
     */
    file("${workflow.workDir}/GEMmaker/process/DONE").delete()

    /**
     * Move any incomplete samples from previous run back to staging.
     */
    file("${workflow.workDir}/GEMmaker/process/*").each { sample_file ->
        sample_file.moveTo("${workflow.workDir}/GEMmaker/stage")
    }

    /**
     * Get list of samples to skip.
     */
    skip_samples = params.skip_samples
        ? file(params.skip_samples).readLines().collect { line -> line.trim() }
        : []

    /**
     * Remove samples in the skip list from staging.
     */
    skip_samples.each { sample_id ->
        skip_sample = file("${workflow.workDir}/GEMmaker/stage/${sample_id}.sample.csv")
        if (skip_sample.exists()) {
            skip_sample.delete()
        }
    }

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
        retrieve_sra_metadata()

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

        hisat2_fpkm_tpm(STRINGTIE_FILES)
        RAW_FILES = hisat2_fpkm_tpm.out.RAW_FILES
        TPM_FILES = hisat2_fpkm_tpm.out.TPM_FILES
        FPKM_FILES = hisat2_fpkm_tpm.out.FPKM_FILES

        // prepare channels for downstream processes
        COMPLETED_SAMPLES = hisat2_fpkm_tpm.out.DONE_SIGNAL

        MULTIQC_FILES = Channel.empty()
            .mix(
                fastqc_1.out.REPORTS.flatten(),
                trimmomatic.out.LOGS.map { it[1] },
                fastqc_2.out.REPORTS.flatten(),
                hisat2.out.LOGS.map { it[1] },
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
     * Clean trimmed fastq files after they are used by fastqc_2 and hisat2.
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
     * Clean sam files after they are used by samtools_sort.
     */
    if ( hisat2_enable && !params.hisat2_keep_sam ) {
        CLEAN_SAM_FILES = SAM_FILES
            .mix(samtools_sort.out.DONE_SIGNAL)
            .groupTuple(size: 2)
            .map { [it[0], it[1].flatten().findAll { v -> v != DONE_SENTINEL }] }

        CLEAN_WORK_FILES = CLEAN_WORK_FILES.mix(CLEAN_SAM_FILES)
    }

    /**
     * Clean bam files after they are used by stringtie.
     */
    if ( hisat2_enable && !params.hisat2_keep_bam ) {
        CLEAN_BAM_FILES = BAM_FILES
            .mix(stringtie.out.DONE_SIGNAL)
            .groupTuple(size: 2)
            .map { [it[0], it[1].flatten().findAll { v -> v != DONE_SENTINEL }] }

        CLEAN_WORK_FILES = CLEAN_WORK_FILES.mix(CLEAN_BAM_FILES)
    }

    /**
     * Clean stringtie files after they are used by hisat2_fpkm_tpm.
     */
    if ( hisat2_enable && !params.hisat2_keep_data ) {
        CLEAN_STRINGTIE_FILES = STRINGTIE_FILES
            .mix(hisat2_fpkm_tpm.out.DONE_SIGNAL)
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

/**
 * Retrieves metadata for all of the remote samples
 * and maps SRA runs to SRA experiments.
 */
process retrieve_sra_metadata {
  publishDir params.outdir, mode: params.publish_dir_mode, pattern: "failed_runs.metadata.txt"
  container "systemsgenetics/gemmaker:2.0.0"

  output:
  stdout emit: SRR2SRX
  path("failed_runs.metadata.txt"), emit: FAILED_RUNS

  script:
  """
  >&2 echo "#TRACE n_remote_run_ids=`cat ${file(params.sras)} | wc -l`"

  retrieve_sra_metadata.py \
      --run_id_file ${file(params.sras)} \
      --meta_dir ${workflow.workDir}/GEMmaker \
       ${params.skip_samples ? "--skip_file ${file(params.skip_samples)}" : ""}
  """
}



/**
 * Move a new sample into the process directory when
 * a previous sample is completed.
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
    sample_file = file("${workflow.workDir}/GEMmaker/process/${sample_id}.sample.csv")
    sample_file.moveTo("${workflow.workDir}/GEMmaker/done")

    // Move the next sample file into the processing directory.
    staged_files = file("${workflow.workDir}/GEMmaker/stage/*")
    if (staged_files.size() > 0) {
        staged_files.first().moveTo("${workflow.workDir}/GEMmaker/process")
    }

    // Write the "done" file if there are no more samples to process.
    else {
        done_file = file("${workflow.workDir}/GEMmaker/process/DONE")
        done_file << ""
    }
}



/**
 * Downloads SRA files from NCBI using the SRA Toolkit.
 */
process download_runs {
  tag { sample_id }
  publishDir params.outdir, mode: params.publish_dir_mode, pattern: '*.failed_runs.download.txt', saveAs: { "Samples/${sample_id}/${it}" }
  container "systemsgenetics/gemmaker:2.0.0"

  input:
    tuple val(sample_id), val(run_ids), val(type)

  output:
    tuple val(sample_id), path("*.sra"), optional: true, emit: SRA_FILES
    tuple val(sample_id), path("sample_failed"), optional: true, emit: FAILED_SAMPLES
    tuple val(sample_id), path("*.failed_runs.download.txt"), emit: FAILED_RUNS

  script:
  """
  echo "#TRACE sample_id=${sample_id}"
  echo "#TRACE n_remote_run_ids=${run_ids.tokenize(" ").size()}"
  echo "#TRACE n_spots=`retrieve_sra_spots.py ${workflow.workDir}/GEMmaker ${sample_id}`"

  retrieve_sra.py --sample ${sample_id} --run_ids ${run_ids} --akey \${ASPERA_KEY}
  """
}



/**
 * Extracts FASTQ files from downloaded SRA files.
 */
process fastq_dump {
  tag { sample_id }
  publishDir params.outdir, mode: params.publish_dir_mode, pattern: publish_pattern_fastq_dump, saveAs: { "Samples/${sample_id}/${it}" }
  publishDir params.outdir, mode: params.publish_dir_mode, pattern: '*.failed_runs.fastq-dump.txt', saveAs: { "Samples/${sample_id}/${it}" }
  container "systemsgenetics/gemmaker:2.0.0"

  input:
    tuple val(sample_id), path(sra_files)

  output:
    tuple val(sample_id), path("*.fastq"), optional: true, emit: FASTQ_FILES
    tuple val(sample_id), path("sample_failed"), optional: true, emit: FAILED_SAMPLES
    tuple val(sample_id), path("*.failed_runs.fastq-dump.txt"), emit: FAILED_RUNS
    tuple val(sample_id), val(DONE_SENTINEL), emit: DONE_SIGNAL

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
  container "systemsgenetics/gemmaker:2.0.0"

  input:
    tuple val(sample_id), path(fastq_files)

  output:
    tuple val(sample_id), path("${sample_id}_?.fastq"), emit: FASTQ_FILES
    tuple val(sample_id), val(DONE_SENTINEL), emit: DONE_SIGNAL

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
  tag { sample_id }
  publishDir "${params.outdir}/Samples/${sample_id}", mode: params.publish_dir_mode, pattern: "*_fastqc.*"
  container "systemsgenetics/gemmaker:2.0.0"

  input:
    tuple val(sample_id), path(fastq_files)

  output:
    tuple path("*_fastqc.html"), path("*_fastqc.zip"), emit: REPORTS
    tuple val(sample_id), val(DONE_SENTINEL), emit: DONE_SIGNAL

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
  tag { sample_id }
  publishDir "${params.outdir}/Samples/${sample_id}", mode: params.publish_dir_mode, pattern: publish_pattern_kallisto_ga
  container "systemsgenetics/gemmaker:2.0.0"

  input:
    tuple val(sample_id), path(fastq_files)
    path(kallisto_index)

  output:
    tuple val(sample_id), path("*.ga", type: "dir"), emit: GA_FILES
    tuple val(sample_id), path("*.kallisto.log"), emit: LOGS
    tuple val(sample_id), val(DONE_SENTINEL), emit: DONE_SIGNAL

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
  tag { sample_id }
  publishDir "${params.outdir}/Samples/${sample_id}", mode: params.publish_dir_mode
  container "systemsgenetics/gemmaker:2.0.0"

  input:
    tuple val(sample_id), path(ga_file)

  output:
    path("*.tpm"), optional: true, emit: TPM_FILES
    path("*.raw"), optional: true, emit: RAW_FILES
    tuple val(sample_id), val(DONE_SENTINEL), emit: DONE_SIGNAL

  script:
  """
  echo "#TRACE sample_id=${sample_id}"

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
  tag { sample_id }
  publishDir "${params.outdir}/Samples/${sample_id}", mode: params.publish_dir_mode, pattern: publish_pattern_salmon_ga
  container "systemsgenetics/gemmaker:2.0.0"

  input:
    tuple val(sample_id), path(fastq_files)
    path(salmon_index)

  output:
    tuple val(sample_id), path("*.ga", type: "dir"), emit: GA_FILES
    tuple val(sample_id), path("*.ga", type: "dir"), emit: LOGS
    tuple val(sample_id), val(DONE_SENTINEL), emit: DONE_SIGNAL

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
  tag { sample_id }
  publishDir "${params.outdir}/Samples/${sample_id}", mode: params.publish_dir_mode
  container "systemsgenetics/gemmaker:2.0.0"

  input:
    tuple val(sample_id), path(ga_file)

  output:
    path("*.Salmon.tpm"), optional: true, emit: TPM_FILES
    path("*.Salmon.raw"), optional: true, emit: RAW_FILES
    tuple val(sample_id), val(DONE_SENTINEL), emit: DONE_SIGNAL

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
  tag { sample_id }
  label "process_medium"
  publishDir "${params.outdir}/Samples/${sample_id}", mode: params.publish_dir_mode, pattern: publish_pattern_trimmomatic
  container "systemsgenetics/gemmaker:2.0.0"

  input:
    tuple val(sample_id), path(fastq_files)
    path(fasta_adapter)

  output:
    tuple val(sample_id), path("*_trim.fastq"), emit: FASTQ_FILES
    tuple val(sample_id), path("*.trim.log"), emit: LOGS
    tuple val(sample_id), val(DONE_SENTINEL), emit: DONE_SIGNAL

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
  tag { sample_id }
  publishDir "${params.outdir}/Samples/${sample_id}", mode: params.publish_dir_mode, pattern: "*_fastqc.*"
  container "systemsgenetics/gemmaker:2.0.0"

  input:
    tuple val(sample_id), path(fastq_files)

  output:
    tuple path("*_fastqc.html"), path("*_fastqc.zip"), emit: REPORTS
    tuple val(sample_id), val(DONE_SENTINEL), emit: DONE_SIGNAL

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
  tag { sample_id }
  label "process_medium"
  publishDir "${params.outdir}/Samples/${sample_id}", mode: params.publish_dir_mode, pattern: "*.log"
  container "systemsgenetics/gemmaker:2.0.0"

  input:
    tuple val(sample_id), path(fastq_files)
    path(indexes)

  output:
    tuple val(sample_id), path("*.sam"), emit: SAM_FILES
    tuple val(sample_id), path("*.sam.log"), emit: LOGS
    tuple val(sample_id), val(DONE_SENTINEL), emit: DONE_SIGNAL

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
  tag { sample_id }
  publishDir "${params.outdir}/Samples/${sample_id}", mode: params.publish_dir_mode, pattern: publish_pattern_samtools_sort
  container "systemsgenetics/gemmaker:2.0.0"

  input:
    tuple val(sample_id), path(sam_file)

  output:
    tuple val(sample_id), path("*.bam"), emit: BAM_FILES
    tuple val(sample_id), val(DONE_SENTINEL), emit: DONE_SIGNAL

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
  tag { sample_id }
  publishDir "${params.outdir}/Samples/${sample_id}", mode: params.publish_dir_mode, pattern: publish_pattern_samtools_index
  container "systemsgenetics/gemmaker:2.0.0"

  input:
    tuple val(sample_id), path(bam_file)

  output:
    tuple val(sample_id), path(bam_file), emit: BAM_FILES
    tuple val(sample_id), path("*.bam.bai"), emit: BAI_FILES
    tuple val(sample_id), path("*.bam.log"), emit: LOGS

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
  tag { sample_id }
  label "process_medium"
  publishDir "${params.outdir}/Samples/${sample_id}", mode: params.publish_dir_mode, pattern: publish_pattern_stringtie_ga_gtf
  container "systemsgenetics/gemmaker:2.0.0"

  input:
    tuple val(sample_id), path(bam_file)
    path(gtf_file)

  output:
    tuple val(sample_id), path("*.Hisat2.*"), emit: GA_GTF_FILES
    tuple val(sample_id), val(DONE_SENTINEL), emit: DONE_SIGNAL

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
  tag { sample_id }
  publishDir "${params.outdir}/Samples/${sample_id}", mode: params.publish_dir_mode
  container "systemsgenetics/gemmaker:2.0.0"

  input:
    tuple val(sample_id), path(stringtie_files)

  output:
    path("*.Hisat2.fpkm"), optional: true, emit: FPKM_FILES
    path("*.Hisat2.tpm"), optional: true, emit: TPM_FILES
    path("*.Hisat2.raw"), optional: true, emit: RAW_FILES
    tuple val(sample_id), val(DONE_SENTINEL), emit: DONE_SIGNAL

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
  cache false
  publishDir "${params.outdir}/reports", mode: params.publish_dir_mode
  container "systemsgenetics/gemmaker:2.0.0"

  input:
    val(signal)
    path(config_file)
    path(custom_logo)
    path(input_files)

  output:
    path("multiqc_data"), emit: MULTIQC_DATA
    path("multiqc_report.html"), emit: MULTIQC_REPORT

  script:
  """
  multiqc --config ${config_file} ./
  """
}



/**
 * Creates the GEM file from all the FPKM/TPM outputs
 */
process create_gem {
  publishDir "${params.outdir}/GEMs", mode: params.publish_dir_mode
  container "systemsgenetics/gemmaker:2.0.0"

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
  publishDir "${params.outdir}/reports", mode: params.publish_dir_mode
  container "systemsgenetics/gemmaker:2.0.0"

  input:
    path(failed_runs)
    path(failed_run_template)

  output:
    path("failed_SRA_run_report.html"), emit: FAILED_RUN_REPORT

  script:
  """
  failed_runs_report.py --template ${failed_run_template}
  """
}



/**
 * Clean intermediate files.
 */
process clean_work_files {
  tag { sample_id }
  label "local"

  input:
    tuple val(sample_id), val(files)

  script:
  """
  clean_work_files.sh "${files.join(" ")}"
  """
}



/**
 * Clean intermediate directories.
 */
process clean_work_dirs {
  tag { sample_id }
  label "local"

  input:
    tuple val(sample_id), val(directory)

  script:
  """
  clean_work_dirs.sh "${directory}"
  """
}
