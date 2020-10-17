#!/usr/bin/env nextflow

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



/**
 * Check to make sure that required reference files exist
 */
// If Hisat2 was selected:
if (selected_tool == 0) {
  gtfFile = file("${params.input.reference_dir}/${params.input.hisat2.gtf_file}")

  if (gtfFile.isEmpty()) {
    error "Error: GTF reference file for Hisat2 does not exist or is empty! Please Check that you have the proper references, that they are placed in the reference directory, and they are named properly.\
    \nGEMmaker is missing the following file: '${params.input.reference_dir}/${params.input.hisat2.gtf_file}' (where '*' is the name of your organism)"
  }

  hisat2_index_dir = file("${params.input.reference_dir}/${params.input.hisat2.index_dir}")

  if (!hisat2_index_dir.isDirectory()) {
    error "Error: hisat2 Index Directory does not exist or is empty! Please Check that you have the proper references, that they are placed in the reference directory, and they are named properly.\
    \nGEMmaker is missing the following file: '${params.input.reference_dir}/${params.input.hisat2.index_dir}' (where '*' is the name of your organism)"
  }
}

// If Kallisto was selected
if (selected_tool == 1) {
  kallisto_index_file = file("${params.input.reference_dir}/${params.input.kallisto.index_file}")

  if (kallisto_index_file.isEmpty()) {
    error "Error: Kallisto Index File does not exist or is empty! Please Check that you have the proper references, that they are placed in the reference directory, and they are named properly.\
    \nGEMmaker is missing the following file: '${params.input.reference_dir}/${params.input.kallisto.index_file}' (where '*' is the name of your organism)"
  }
}

// If Salmon was selected
if (selected_tool == 2) {
  salmon_index_dir = file("${params.input.reference_dir}/${params.input.salmon.index_dir}")

  if (!salmon_index_dir.isDirectory()) {
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
 * Local Sample Input.
 * This checks the folder that the user has given
 */
if (params.input.local_samples_path == "none") {
  LOCAL_SAMPLE_FILES = Channel.empty()
}
else {
  LOCAL_SAMPLE_FILES = Channel.fromFilePairs( "${params.input.input_data_dir}/${params.input.local_sample_files}", size: -1 )
}

/**
 * Remote fastq_run_id Input.
 */
if (params.input.remote_list_path == "none") {
  SRR_FILE = Channel.empty()
}
else {
  SRR_FILE = Channel.fromPath("${params.input.input_data_dir}/${params.input.remote_sample_list}")
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
  >&2 echo "#TRACE n_remote_run_ids=`cat ${srr_file} | wc -l`"

  retrieve_sra_metadata.py ${srr_file}
  """
}



/**
 * Merge remote samples and local samples into one channel.
 */
REMOTE_SAMPLES_LIST
  .splitCsv()
  .groupTuple(by: 1)
  .map{ [it[1], "remote", it[0], []] }
  .set{REMOTE_SAMPLES}

LOCAL_SAMPLE_FILES
  .map{ [it[0], "local", [], it[1]] }
  .set{LOCAL_SAMPLES}

ALL_SAMPLES = REMOTE_SAMPLES.mix(LOCAL_SAMPLES)

/**
 * Create channel for index files based on the selected aligner.
 */
if ( params.input.hisat2.enable == true ) {
  INDEXES = HISAT2_INDEXES
}
else if ( params.input.kallisto.enable == true ) {
  INDEXES = KALLISTO_INDEX
}
else if ( params.input.salmon.enable == true ) {
  INDEXES = SALMON_INDEXES
}



/**
 * Process a single sample end-to-end.
 *
 * This process requires that the ILLUMINACLIP_PATH environment
 * variable be set in the trimmomatic module. This indicates
 * the path where the clipping files are stored.
 *
 * MINLEN is calculated using based on percentage of the mean
 * read length. The percenage is determined by the user in the
 * "nextflow.config" file
 */
process process_sample {
  tag { sample_id }
  label "gemmaker"
  label "multithreaded"
  label "retry_ignore"
  publishDir params.output.sample_dir, mode: params.output.publish_mode

  input:
    set val(sample_id), val(sample_type), val(remote_ids), val(local_files) from ALL_SAMPLES
    file fasta_adapter from FASTA_ADAPTER
    file indexes from INDEXES
    file gtf_file from GTF_FILE

  output:
    val(sample_id) into COMPLETED_SAMPLES
    file("*.sra")      optional true into SRA_FILES
    file("*.fastq")    optional true into FASTQ_FILES
    file("*_fastqc.*") optional true into FASTQC_FILES
    file("*.log")      optional true into LOG_FILES
    file("*.sam")      optional true into SAM_FILES
    file("*.bam")      optional true into BAM_FILES
    file("*.bam.bai")  optional true into BAI_FILES
    file("*.ga")       optional true into GA_FILES
    file("*.gtf")      optional true into GTF_FILES
    file("*.raw")      optional true into RAW_FILES
    file("*.fpkm")     optional true into FPKM_FILES
    file("*.tpm")      optional true into TPM_FILES

  script:
  """
  echo "#TRACE sample_id=${sample_id}"
  echo "#TRACE sample_type=${sample_type}"
  echo "#TRACE hisat2_enable=${params.input.hisat2.enable}}"
  echo "#TRACE kallisto_enable=${params.input.kallisto.enable}}"
  echo "#TRACE salmon_enable=${params.input.salmon.enable}}"
  echo "#TRACE publish_sra=${params.output.publish_sra}}"
  echo "#TRACE publish_downloaded_fastq=${params.output.publish_downloaded_fastq}}"
  echo "#TRACE publish_trimmed_fastq=${params.output.publish_trimmed_fastq}}"
  echo "#TRACE publish_bam=${params.output.publish_bam}}"
  echo "#TRACE publish_stringtie_gtf_and_ga=${params.output.publish_stringtie_gtf_and_ga}}"
  echo "#TRACE publish_gene_abundance=${params.output.publish_gene_abundance}}"
  echo "#TRACE publish_fpkm=${params.output.publish_fpkm}}"
  echo "#TRACE publish_tpm=${params.output.publish_tpm}}"
  echo "#TRACE publish_raw=${params.output.publish_raw}}"

  # for remote samples, prepare FASTQ files from NCBI
  if [[ "${sample_type}" == "remote" ]]; then
    # download SRA files from NCBI
    echo "#TRACE n_remote_run_ids=${remote_ids.size()}"

    retrieve_sra.py ${remote_ids.join(',')}

    # extract FASTQ files from SRA files
    echo "#TRACE sra_bytes=`stat -c '%s' *.sra | paste -sd+ | bc`"

    sra2fastq.py *.sra

    # remove SRA files if they will not be published
    if [[ ${params.output.publish_sra} == false ]]; then
      rm -f *.sra
    fi

    # merge the FASTQ files from each run in the experiment
    DOWNLOADED_FASTQ_FILES=\$(ls *.fastq)

    fastq_merge.sh ${sample_id}

    # remove downloaded FASTQ files if they will not be published
    if [[ ${params.output.publish_downloaded_fastq} == false ]]; then
      rm -f \${DOWNLOADED_FASTQ_FILES}
    fi

  # for local samples, fetch FASTQ files from filesystem
  elif [[ "${sample_type}" == "local" ]]; then
    cp ${local_files.join(' ')} .
  fi

  # perform fastqc on raw FASTQ files
  MERGED_FASTQ_FILES=\$(ls ${sample_id}_?.fastq)
  echo "#TRACE fastq_lines=`cat \${MERGED_FASTQ_FILES} | wc -l`"

  fastqc \${MERGED_FASTQ_FILES}

  # use hisat2 for alignment if specified
  if [[ ${params.input.hisat2.enable} == "true" ]]; then
    # perform trimmomatic on all fastq files
    echo "#TRACE n_cpus=${task.cpus}"
    echo "#TRACE minlen=${params.software.trimmomatic.MINLEN}"
    echo "#TRACE leading=${params.software.trimmomatic.LEADNING}"
    echo "#TRACE trailing=${params.software.trimmomatic.TRAILING}"
    echo "#TRACE slidingwindow=${params.software.trimmomatic.SLIDINGWINDOW}"
    echo "#TRACE fasta_lines=`cat ${fasta_adapter} | wc -l`"

    trimmomatic.sh \
      ${sample_id} \
      ${params.software.trimmomatic.MINLEN} \
      ${task.cpus} \
      ${fasta_adapter} \
      ${params.software.trimmomatic.LEADING} \
      ${params.software.trimmomatic.TRAILING} \
      ${params.software.trimmomatic.SLIDINGWINDOW} \
      ${params.software.trimmomatic.quality}

    # remove merged fastq files if they will not be published
    if [[ ${params.output.publish_downloaded_fastq} == false ]]; then
      rm -f \${MERGED_FASTQ_FILES}
    fi

    # perform fastqc on all trimmed fastq files
    TRIMMED_FASTQ_FILES=\$(ls ${sample_id}_*trim.fastq)
    echo "#TRACE trimmed_fastq_lines=`cat \${TRIMMED_FASTQ_FILES} | wc -l`"

    fastqc \${TRIMMED_FASTQ_FILES}

    # perform hisat2 alignment of fastq files to a genome reference
    echo "#TRACE index_bytes=`stat -c '%s' ${indexes} | paste -sd+ | bc`"

    hisat2.sh \
      ${sample_id} \
      ${params.input.reference_name} \
      ${task.cpus}

    # remove unused hisat2 output files
    rm -f ${sample_id}_un.fastq

    # remove trimmed fastq files if they will not be published
    if [[ ${params.output.publish_trimmed_fastq} == false ]]; then
      rm -f \${TRIMMED_FASTQ_FILES}
    fi

    # sort the SAM alignment file and convert it to BAM
    echo "#TRACE sam_lines=`cat *.sam | wc -l`"

    samtools sort \
      -o ${sample_id}_vs_${params.input.reference_name}.bam \
      -O bam \
      -T temp \
      ${sample_id}_vs_${params.input.reference_name}.sam

    # remove SAM file as it will not be published
    rm -f *.sam

    # index BAM alignment file
    echo "#TRACE bam_bytes=`stat -c '%s' *.bam`"

    samtools index ${sample_id}_vs_${params.input.reference_name}.bam
    samtools stats ${sample_id}_vs_${params.input.reference_name}.bam > ${sample_id}_vs_${params.input.reference_name}.bam.log

    # generate expression-level transcript abundance
    echo "#TRACE gtf_lines=`cat *.gtf | wc -l`"

    stringtie \
      -v \
      -p ${task.cpus} \
      -e \
      -o ${sample_id}_vs_${params.input.reference_name}.Hisat2.gtf \
      -G ${gtf_file} \
      -A ${sample_id}_vs_${params.input.reference_name}.Hisat2.ga \
      -l ${sample_id} \
      ${sample_id}_vs_${params.input.reference_name}.bam

    # remove BAM file if it will not be published
    if [[ ${params.output.publish_bam} == false ]]; then
      rm -f *.bam
      rm -f *.bam.bai
    fi

    # generate raw counts from hisat2/stringtie
    echo "#TRACE ga_lines=`cat *.ga | wc -l`"
    echo "#TRACE gtf_lines=`cat *.gtf | wc -l`"

    hisat2_fpkm_tpm.sh \
      ${params.output.publish_fpkm} \
      ${sample_id} \
      ${params.input.reference_name} \
      ${params.output.publish_tpm} \
      ${params.output.publish_raw}

    # remove stringtie output files if they will not be published
    if [[ ${params.output.publish_stringtie_gtf_and_ga} == false ]]; then
      rm -rf *.ga
      rm -rf *.gtf
    fi

  # or use KALLISTO if specified
  elif [[ ${params.input.kallisto.enable} == "true" ]]; then
    # perform KALLISTO alignment of fastq files
    echo "#TRACE index_bytes=`stat -c '%s' ${indexes}`"

    kallisto.sh \
      ${sample_id} \
      ${indexes} \
      ${params.input.reference_name}

    # remove merged fastq files if they will not be published
    if [[ ${params.output.publish_downloaded_fastq} == false ]]; then
      rm -f \${MERGED_FASTQ_FILES}
    fi

    # generate final TPM and raw count files
    echo "#TRACE ga_lines=`cat *.Kallisto.ga | wc -l`"

    kallisto_tpm.sh \
      ${sample_id} \
      ${params.output.publish_tpm} \
      ${params.output.publish_raw} \
      ${params.input.reference_name}

    # remove gene abundance files if they will not be published
    if [[ ${params.output.publish_gene_abundance} == false ]]; then
      rm -rf *.ga
    fi

  # or use SALMON if specified
  elif [[ ${params.input.salmon.enable} == "true" ]]; then
    # perform SALMON alignment of fastq files
    echo "#TRACE index_bytes=`stat -c '%s' ${indexes} | paste -sd+ | bc`"

    salmon.sh \
      ${sample_id} \
      ${task.cpus} \
      ${params.input.reference_name}

    # remove merged fastq files if they will not be published
    if [[ ${params.output.publish_downloaded_fastq} == false ]]; then
      rm -f \${MERGED_FASTQ_FILES}
    fi

    # generate final TPM and raw count files
    salmon_tpm.sh \
      ${params.output.publish_tpm} \
      ${sample_id} \
      ${params.input.reference_name} \
      ${params.output.publish_raw}

    # remove gene abundance files if they will not be published
    if [[ ${params.output.publish_gene_abundance} == false ]]; then
      rm -rf `find *.ga -type f | egrep -v "aux_info/meta_info.json|/libParams/flenDist.txt"`
    fi
  fi
  """
}



/**
 * Send completed samples to each process that uses them
 */
COMPLETED_SAMPLES.into { COMBINED_SAMPLES_FOR_MULTIQC; COMPLETED_SAMPLES_FOR_GEM }



/**
 * Process to generate the multiqc report once everything is completed
 */
process multiqc {
  label "multiqc"
  publishDir "${params.output.dir}/reports", mode: params.output.publish_mode

  input:
    val signal from COMBINED_SAMPLES_FOR_MULTIQC.collect()

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
 * Creates the GEM file from all the FPKM/TPM outputs
 */
process create_gem {
  label "python3"
  publishDir "${params.output.dir}/GEMs", mode: params.output.publish_mode

  input:
    val signal from COMPLETED_SAMPLES_FOR_GEM.collect()

  output:
    file "*.GEM.*.txt" into GEM_FILES

  when:
    params.output.create_gem == true

  script:
  """
  echo "#TRACE publish_fpkm=${params.output.publish_fpkm}}"
  echo "#TRACE hisat2_enable=${params.input.hisat2.enable}}"
  echo "#TRACE publish_tpm=${params.output.publish_tpm}}"
  echo "#TRACE publish_raw=${params.output.publish_raw}}"
  echo "#TRACE fpkm_lines=`cat ${workflow.launchDir}/${params.output.dir}/*.fpkm 2> /dev/null | wc -l`"
  echo "#TRACE tpm_lines=`cat ${workflow.launchDir}/${params.output.dir}/*.tpm 2> /dev/null | wc -l`"
  echo "#TRACE raw_lines=`cat ${workflow.launchDir}/${params.output.dir}/*.raw 2> /dev/null | wc -l`"

  create_gem.sh \
    ${params.output.publish_fpkm} \
    ${params.input.hisat2.enable} \
    ${workflow.launchDir} \
    ${params.output.dir} \
    ${params.project.machine_name} \
    ${params.output.publish_raw} \
    ${params.output.publish_tpm}
  """
}
