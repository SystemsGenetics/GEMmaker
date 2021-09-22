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
    method_lock_file << params.pipeline
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



workflow {
    get_software_versions()

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

    LOCAL_SAMPLES = LOCAL_SAMPLES
        .map{ [it[0], it[1], "local"] }

    /**
     * Retrieve metadata for remote samples from NCBI SRA.
     */
    SRR_FILE = (params.sras != "")
        ? Channel.fromPath( params.sras )
        : Channel.empty()

    SKIP_SAMPLES_FILE = params.skip_samples
        ? Channel.fromPath( params.skip_samples )
        : Channel.fromList( [null] )

    retrieve_sra_metadata(SRR_FILE, SKIP_SAMPLES_FILE)

    /**
     * Parse remote samples from the SRR2SRX mapping.
     */
    REMOTE_SAMPLES = retrieve_sra_metadata.out.SRR2SRX
        .splitCsv()
        .groupTuple(by: 1)
        .map { [it[1], it[0].join(" "), "remote"] }

    /**
     * Prepare sample files for staging and processing.
     */
    LOCAL_SAMPLES.mix(REMOTE_SAMPLES)
        .filter { sample -> !skip_samples.contains(sample[0]) }
        .subscribe onNext: { sample ->
            // Stage samples that are not in the skip list.
            (sample_id, run_files_or_ids, sample_type) = sample

            sample_file = file("${workflow.workDir}/GEMmaker/stage/${sample_id}.sample.csv")

            line = (sample_type == "local")
                ? "${sample_id}\t${run_files_or_ids}\t${sample_type}"
                : "${sample_id}\t${run_files_or_ids}\t${sample_type}"
            sample_file << line
        }, onComplete: {
            // Move the first batch of samples into processing.
            file("${workflow.workDir}/GEMmaker/stage/*.sample.csv")
                .sort()
                .take(params.max_cpus)
                .each { sample ->
                    sample.moveTo("${workflow.workDir}/GEMmaker/process")
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
        .join(LOCAL_SAMPLES)
        .map { [it[0], it[3]] }

    /**
     * Download, extract, and merge remote samples.
     */
    REMOTE_SAMPLES = NEXT_SAMPLE.remote

    download_runs(REMOTE_SAMPLES)
    SRA_FILES = download_runs.out.SRA_FILES

    fastq_dump(SRA_FILES)
    DOWNLOADED_FASTQ_FILES = fastq_dump.out.FASTQ_FILES

    fastq_merge(DOWNLOADED_FASTQ_FILES)
    MERGED_FASTQ_FILES = fastq_merge.out.FASTQ_FILES

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
            download_runs.out.FAILED_SAMPLES.map { it[0] },
            fastq_dump.out.FAILED_SAMPLES.map { it[0] })

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
     * completed or if on a resume when the bootstrap signal is
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
    FAILED_RUNS = Channel.empty()
        .mix(
            retrieve_sra_metadata.out.FAILED_RUNS,
            download_runs.out.FAILED_RUNS.map { it[1] },
            fastq_dump.out.FAILED_RUNS.map { it[1] })
    FAILED_RUN_TEMPLATE = Channel.fromPath( params.failed_run_report_template )

    failed_run_report(
        FAILED_RUNS.collect(),
        FAILED_RUN_TEMPLATE.collect())

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

    /**
     * Clean sra files after they are used by fastq_dump.
     */
    if ( !params.keep_sra ) {
        CLEAN_SRA_READY = SRA_FILES
            .mix(fastq_dump.out.DONE_SIGNAL)
            .groupTuple(size: 2)
            .map { [it[0], it[1].flatten().findAll { v -> v != DONE_SENTINEL }] }

        clean_sra(CLEAN_SRA_READY)
    }

    /**
     * Clean downloaded fastq files after they are used by fastq_merge.
     */
    if ( !params.keep_retrieved_fastq ) {
        CLEAN_DOWNLOADED_FASTQ_READY = DOWNLOADED_FASTQ_FILES
            .mix(fastq_merge.out.DONE_SIGNAL)
            .groupTuple(size: 2)
            .map { [it[0], it[1].flatten().findAll { v -> v != DONE_SENTINEL }] }

        clean_downloaded_fastq(CLEAN_DOWNLOADED_FASTQ_READY)
    }

    /**
     * Clean merged fastq files after they are used by fastqc_1 and
     * the quantification tool.
     */
    if ( !params.keep_retrieved_fastq ) {
        CLEAN_MERGED_FASTQ_READY = MERGED_FASTQ_FILES
            .mix(
                fastqc_1.out.DONE_SIGNAL,
                MERGED_FASTQ_DONE)
            .groupTuple(size: 3)
            .map { [it[0], it[1].flatten().findAll { v -> v != DONE_SENTINEL }] }

        clean_merged_fastq(CLEAN_MERGED_FASTQ_READY)
    }

    /**
     * Clean trimmed fastq files after they are used by fastqc_2 and hisat2.
     */
    if ( hisat2_enable && !params.trimmomatic_keep_trimmed_fastq ) {
        CLEAN_TRIMMED_FASTQ_READY = TRIMMED_FASTQ_FILES
            .mix(
                fastqc_2.out.DONE_SIGNAL,
                hisat2.out.DONE_SIGNAL)
            .groupTuple(size: 3)
            .map { [it[0], it[1].flatten().findAll { v -> v != DONE_SENTINEL }] }

        clean_trimmed_fastq(CLEAN_TRIMMED_FASTQ_READY)
    }

    /**
     * Clean sam files after they are used by samtools_sort.
     */
    if ( hisat2_enable && !params.hisat2_keep_sam ) {
        CLEAN_SAM_READY = SAM_FILES
            .mix(samtools_sort.out.DONE_SIGNAL)
            .groupTuple(size: 2)
            .map { [it[0], it[1].flatten().findAll { v -> v != DONE_SENTINEL }] }

        clean_sam(CLEAN_SAM_READY)
    }

    /**
     * Clean bam files after they are used by stringtie.
     */
    if ( hisat2_enable && !params.hisat2_keep_bam ) {
        CLEAN_BAM_READY = BAM_FILES
            .mix(stringtie.out.DONE_SIGNAL)
            .groupTuple(size: 2)
            .map { [it[0], it[1].flatten().findAll { v -> v != DONE_SENTINEL }] }

        clean_bam(CLEAN_BAM_READY)
    }

    /**
     * Clean stringtie files after they are used by hisat2_fpkm_tpm.
     */
    if ( hisat2_enable && !params.hisat2_keep_data ) {
        CLEAN_STRINGTIE_READY = STRINGTIE_FILES
            .mix(hisat2_fpkm_tpm.out.DONE_SIGNAL)
            .groupTuple(size: 2)
            .map { [it[0], it[1].flatten().findAll { v -> v != DONE_SENTINEL }] }

        clean_stringtie_ga(CLEAN_STRINGTIE_READY)
    }

    /**
     * Clean kallisto files after they are used by kallisto_tpm.
     */
    if ( kallisto_enable && !params.kallisto_keep_data ) {
        CLEAN_KALLISTO_GA_READY = GA_FILES
            .mix(kallisto_tpm.out.DONE_SIGNAL)
            .groupTuple(size: 2)
            .map { [it[0], it[1][0]] }

        clean_kallisto_ga(CLEAN_KALLISTO_GA_READY)
    }

    /**
     * Clean salmon files after they are used by salmon_tpm.
     */
    if ( salmon_enable && !params.salmon_keep_data ) {
        CLEAN_SALMON_GA_READY = GA_FILES
            .mix(salmon_tpm.out.DONE_SIGNAL)
            .groupTuple(size: 2)
            .map { [it[0], it[1].flatten().findAll { v -> v != DONE_SENTINEL }] }

        clean_salmon_ga(CLEAN_SALMON_GA_READY)
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
    path(skip_file)

  output:
    stdout emit: SRR2SRX
    path("failed_runs.metadata.txt"), emit: FAILED_RUNS

  script:
  """
  >&2 echo "#TRACE n_remote_run_ids=`cat ${srr_file} | wc -l`"

  retrieve_sra_metadata.py \
      --run_id_file ${srr_file} \
      --meta_dir ${workflow.workDir}/GEMmaker \
      ${skip_file != null ? "--skip_file ${skip_file}" : ""}
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
  publishDir params.outdir, mode: params.publish_dir_mode, pattern: '*.failed_runs.download.txt', saveAs: { "Samples/${sample_id}/${it}" }

  tag { sample_id }
  label "download_runs"

  input:
    tuple val(sample_id), val(run_ids), val(type)

  output:
    tuple val(sample_id), path("*.sra"), optional: true, emit: SRA_FILES
    tuple val(sample_id), path("sample_failed"), optional: true, emit: FAILED_SAMPLES
    tuple val(sample_id), path("*.failed_runs.download.txt"), emit: FAILED_RUNS

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
  publishDir "${params.outdir}/Samples/${sample_id}", mode: params.publish_dir_mode, pattern: "*_fastqc.*"
  tag { sample_id }
  label "fastqc"

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
  publishDir "${params.outdir}/Samples/${sample_id}", mode: params.publish_dir_mode, pattern: publish_pattern_kallisto_ga
  tag { sample_id }
  label "process_medium"
  label "kallisto"

  input:
    tuple val(sample_id), path(fastq_files)
    path(kallisto_index)

  output:
    tuple val(sample_id), path("*.ga"), emit: GA_FILES
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
  publishDir "${params.outdir}/Samples/${sample_id}", mode: params.publish_dir_mode
  tag { sample_id }

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
  publishDir "${params.outdir}/Samples/${sample_id}", mode: params.publish_dir_mode, pattern: publish_pattern_salmon_ga
  tag { sample_id }
  label "process_medium"
  label "salmon"

  input:
    tuple val(sample_id), path(fastq_files)
    path(salmon_index)

  output:
    tuple val(sample_id), path("*.ga"), emit: GA_FILES
    tuple val(sample_id), path("*.ga", type: 'dir'), emit: LOGS
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
  publishDir "${params.outdir}/Samples/${sample_id}", mode: params.publish_dir_mode
  tag { sample_id }

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
  publishDir "${params.outdir}/Samples/${sample_id}", mode: params.publish_dir_mode, pattern: publish_pattern_trimmomatic
  tag { sample_id }
  label "process_medium"
  label "trimmomatic"

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
  publishDir "${params.outdir}/Samples/${sample_id}", mode: params.publish_dir_mode, pattern: "*_fastqc.*"
  tag { sample_id }
  label "fastqc"

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
  publishDir "${params.outdir}/Samples/${sample_id}", mode: params.publish_dir_mode, pattern: "*.log"
  tag { sample_id }
  label "process_medium"
  label "hisat2"

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
  publishDir "${params.outdir}/Samples/${sample_id}", mode: params.publish_dir_mode, pattern: publish_pattern_samtools_sort
  tag { sample_id }
  label "samtools"

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
  publishDir "${params.outdir}/Samples/${sample_id}", mode: params.publish_dir_mode, pattern: publish_pattern_samtools_index
  tag { sample_id }
  label "samtools"

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
  publishDir "${params.outdir}/Samples/${sample_id}", mode: params.publish_dir_mode, pattern: publish_pattern_stringtie_ga_gtf
  tag { sample_id }
  label "process_medium"
  label "stringtie"

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
  publishDir "${params.outdir}/Samples/${sample_id}", mode: params.publish_dir_mode
  tag { sample_id }
  label "stringtie"

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
  label "multiqc"
  cache false
  publishDir "${params.outdir}/reports", mode: params.publish_dir_mode

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
 * Clean up downloaded SRA files
 */
process clean_sra {
  tag { sample_id }

  input:
    tuple val(sample_id), val(files)

  script:
  """
  clean_work_files.sh "${files.join(" ")}"
  """
}



/**
 * Clean up downloaded fastq files
 */
process clean_downloaded_fastq {
  tag { sample_id }

  input:
    tuple val(sample_id), val(files)

  script:
  """
  clean_work_files.sh "${files.join(" ")}"
  """
}



/**
 * Clean up merged fastq files
 */
process clean_merged_fastq {
  tag { sample_id }

  input:
    tuple val(sample_id), val(files)

  script:
  """
  clean_work_files.sh "${files.join(" ")}"
  """
}



/**
 * Clean up trimmed fastq files
 */
process clean_trimmed_fastq {
  tag { sample_id }

  input:
    tuple val(sample_id), val(files)

  script:
  """
  clean_work_files.sh "${files.join(" ")}"
  """
}



/**
 * Clean up SAM files
 */
process clean_sam {
  tag { sample_id }

  input:
    tuple val(sample_id), val(files)

  script:
  """
  clean_work_files.sh "${files.join(" ")}"
  """
}



/**
 * Clean up BAM files
 */
process clean_bam {
  tag { sample_id }

  input:
    tuple val(sample_id), val(files)

  script:
  """
  clean_work_files.sh "${files.join(" ")}"
  """
}



/**
 * Clean up stringtie GA files
 */
process clean_stringtie_ga {
  tag { sample_id }

  input:
    tuple val(sample_id), val(files)

  script:
  """
  clean_work_files.sh "${files.join(" ")}"
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
  clean_work_dirs.sh "${directory}"
  """
}



/**
 * Clean up Salmon GA files
 */
process clean_salmon_ga {
  tag { sample_id }

  input:
    tuple val(sample_id), val(files)

  script:
  """
  clean_work_files.sh "${files.join(" ")}"
  """
}
