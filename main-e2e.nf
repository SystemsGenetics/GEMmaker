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
  Remote fastq list path:     ${params.input.remote_list_path}
  Local sample glob:          ${params.input.local_samples_path}
  Reference genome path:      ${params.input.reference_path}
  Reference genome prefix:    ${params.input.reference_prefix}


Output Parameters:
------------------
  Output directory:           ${workflow.launchDir}/${params.output.dir}
  Publish SRA:                ${params.output.publish_sra}
  Publish downloaded FASTQ:   ${params.output.publish_downloaded_fastq}
  Publish trimmed FASTQ:      ${params.output.publish_trimmed_fastq}
  Publish BAM:                ${params.output.publish_bam}
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



/**
 * Create value channels that can be reused
 */
HISAT2_INDEXES = Channel.fromPath("${params.input.reference_path}/${params.input.reference_prefix}*.ht2*").collect()
KALLISTO_INDEX = Channel.fromPath("${params.input.reference_path}/${params.input.reference_prefix}.transcripts.Kallisto.indexed").collect()
SALMON_INDEXES = Channel.fromPath("${params.input.reference_path}/${params.input.reference_prefix}.transcripts.Salmon.indexed/*").collect()
FASTA_ADAPTER = Channel.fromPath("${params.software.trimmomatic.clip_path}").collect()
GTF_FILE = Channel.fromPath("${params.input.reference_path}/${params.input.reference_prefix}.gtf").collect()



/**
 * Local Sample Input.
 * This checks the folder that the user has given
 */
if (params.input.local_samples_path == "none") {
  LOCAL_SAMPLE_FILES = Channel.empty()
}
else {
  LOCAL_SAMPLE_FILES = Channel.fromFilePairs( params.input.local_samples_path, size: -1 )
}

/**
 * Remote fastq_run_id Input.
 */
if (params.input.remote_list_path == "none") {
  SRR_FILE = Channel.empty()
}
else {
  SRR_FILE = Channel.fromPath(params.input.remote_list_path)
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
if ( params.software.alignment == 0 ) {
  INDEXES = HISAT2_INDEXES
}
else if ( params.software.alignment == 1 ) {
  INDEXES = KALLISTO_INDEX
}
else if ( params.software.alignment == 2 ) {
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
  publishDir params.output.sample_dir, mode: params.output.publish_mode

  input:
    set val(sample_id), val(type), val(remote_ids), val(local_files) from ALL_SAMPLES
    file fasta_adapter from FASTA_ADAPTER
    file indexes from INDEXES
    file gtf_file from GTF_FILE

  output:
    val(sample_id) into COMPLETED_SAMPLES
    file("*.sra") optional true into SRA_FILES
    file("*.fastq") optional true into FASTQ_FILES
    file("*fastqc.*") optional true into FASTQC_FILES
    file("*.log") optional true into LOG_FILES
    file("*.sam") optional true into SAM_FILES
    file("*.bam") optional true into BAM_FILES
    file("*.bam.bai") optional true into BAI_FILES
    file("*.ga") optional true into GA_FILES
    file("*.gtf") optional true into GTF_FILES
    file("*.raw") optional true into RAW_FILES
    file("*.fpkm") optional true into FPKM_FILES
    file("*.tpm") optional true into TPM_FILES

  script:
  """
  # for remote samples, prepare FASTQ files from NCBI
  if [[ "${type}" == "remote" ]]; then
    # download SRA files from NCBI
    SRR_IDS="${remote_ids.join(' ')}"

    for id in \$SRR_IDS; do
      ascp_path=`which ascp`
      prefetch -v --max-size 50G --output-directory . --ascp-path "\$ascp_path|\$ASPERA_KEY" --ascp-options "-k 1 -T -l 1000m" \$id
    done

    # extract FASTQ files from SRA files
    SRA_FILES=\$(ls *.sra)

    for sra_file in \$SRA_FILES; do
      fastq-dump --split-files \$sra_file
    done

    # remove SRA files if they will not be published
    if [[ ${params.output.publish_sra} == false ]]; then
      rm -f \$SRA_FILES
    fi

    # merge the FASTQ files from each run in the experiment
    DOWNLOADED_FASTQ_FILES=\$(ls *.fastq)

    if ls *_1.fastq >/dev/null 2>&1; then
      cat *_1.fastq >> "${sample_id}_1.fastq"
    fi

    if ls *_2.fastq >/dev/null 2>&1; then
      cat *_2.fastq >> "${sample_id}_2.fastq"
    fi

    # remove downloaded FASTQ files if they will not be published
    if [[ ${params.output.publish_downloaded_fastq} == false ]]; then
      rm -f \$DOWNLOADED_FASTQ_FILES
    fi

  # for local samples, fetch FASTQ files from filesystem
  elif [[ "${type}" == "local" ]]; then
    cp ${local_files.join(' ')} .
  fi

  # perform fastqc on raw FASTQ files
  MERGED_FASTQ_FILES=\$(ls ${sample_id}_?.fastq)

  fastqc \$MERGED_FASTQ_FILES

  # use hisat2 for alignment
  if [[ ${params.software.alignment} == 0 ]]; then
    # perform trimmomatic on all fastq files
    # This script calculates average length of fastq files.
    total=0

    # This if statement checks if the data is single or paired data, and checks length accordingly
    # This script returns 1 number, which can be used for the minlen in trimmomatic
    if [ -e ${sample_id}_1.fastq ] && [ -e ${sample_id}_2.fastq ]; then
      for fastq in ${sample_id}_1.fastq ${sample_id}_2.fastq; do
        a=`awk 'NR%4 == 2 {lengths[length(\$0)]++} END {for (l in lengths) {print l, lengths[l]}}' \$fastq \
        | sort \
        | awk '{ print \$0, \$1*\$2}' \
        | awk '{ SUM += \$3 } { SUM2 += \$2 } END { printf("%.0f", SUM / SUM2 * ${params.software.trimmomatic.MINLEN})} '`
      total=(\$a + \$total)
      done
      total=( \$total / 2 )
      minlen=\$total

    elif [ -e ${sample_id}_1.fastq ]; then
      minlen=`awk 'NR%4 == 2 {lengths[length(\$0)]++} END {for (l in lengths) {print l, lengths[l]}}' ${sample_id}_1.fastq \
        | sort \
        | awk '{ print \$0, \$1*\$2}' \
        | awk '{ SUM += \$3 } { SUM2 += \$2 } END { printf("%.0f", SUM / SUM2 * ${params.software.trimmomatic.MINLEN})} '`
    fi

    if [ -e ${sample_id}_1.fastq ] && [ -e ${sample_id}_2.fastq ]; then
      java -Xmx512m org.usadellab.trimmomatic.Trimmomatic \
        PE \
        -threads ${task.cpus} \
        ${params.software.trimmomatic.quality} \
        ${sample_id}_1.fastq \
        ${sample_id}_2.fastq \
        ${sample_id}_1p_trim.fastq \
        ${sample_id}_1u_trim.fastq \
        ${sample_id}_2p_trim.fastq \
        ${sample_id}_2u_trim.fastq \
        ILLUMINACLIP:${params.software.trimmomatic.clip_path}:2:40:15 \
        LEADING:${params.software.trimmomatic.LEADING} \
        TRAILING:${params.software.trimmomatic.TRAILING} \
        SLIDINGWINDOW:${params.software.trimmomatic.SLIDINGWINDOW} \
        MINLEN:"\$minlen" > ${sample_id}.trim.log 2>&1
    else
      # For ease of the next steps, rename the reverse file to the forward.
      # since these are non-paired it really shouldn't matter.
      if [ -e ${sample_id}_2.fastq ]; then
        mv ${sample_id}_2.fastq ${sample_id}_1.fastq
      fi
      # Now run trimmomatic
      java -Xmx512m org.usadellab.trimmomatic.Trimmomatic \
        SE \
        -threads ${task.cpus} \
        ${params.software.trimmomatic.quality} \
        ${sample_id}_1.fastq \
        ${sample_id}_1u_trim.fastq \
        ILLUMINACLIP:${fasta_adapter}:2:40:15 \
        LEADING:${params.software.trimmomatic.LEADING} \
        TRAILING:${params.software.trimmomatic.TRAILING} \
        SLIDINGWINDOW:${params.software.trimmomatic.SLIDINGWINDOW} \
        MINLEN:"\$minlen" > ${sample_id}.trim.log 2>&1
    fi

    # remove merged fastq files if they will not be published
    if [[ ${params.output.publish_downloaded_fastq} == false ]]; then
      rm -f \$MERGED_FASTQ_FILES
    fi

    # perform fastqc on all trimmed fastq files
    TRIMMED_FASTQ_FILES=\$(ls ${sample_id}_*trim.fastq)

    fastqc \$TRIMMED_FASTQ_FILES

    # perform hisat2 alignment of fastq files to a genome reference
    if [ -e ${sample_id}_2p_trim.fastq ]; then
      hisat2 \
        -x ${params.input.reference_prefix} \
        --no-spliced-alignment \
        -q \
        -1 ${sample_id}_1p_trim.fastq \
        -2 ${sample_id}_2p_trim.fastq \
        -U ${sample_id}_1u_trim.fastq,${sample_id}_2u_trim.fastq \
        -S ${sample_id}_vs_${params.input.reference_prefix}.sam \
        -t \
        -p ${task.cpus} \
        --un ${sample_id}_un.fastq \
        --dta-cufflinks \
        --new-summary \
        --summary-file ${sample_id}_vs_${params.input.reference_prefix}.sam.log
    else
      hisat2 \
        -x ${params.input.reference_prefix} \
        --no-spliced-alignment \
        -q \
        -U ${sample_id}_1u_trim.fastq \
        -S ${sample_id}_vs_${params.input.reference_prefix}.sam \
        -t \
        -p ${task.cpus} \
        --un ${sample_id}_un.fastq \
        --dta-cufflinks \
        --new-summary \
        --summary-file ${sample_id}_vs_${params.input.reference_prefix}.sam.log
    fi

    rm -f ${sample_id}_un.fastq

    # remove trimmed fastq files if they will not be published
    if [[ ${params.output.publish_trimmed_fastq} == false ]]; then
      rm -f \$TRIMMED_FASTQ_FILES
    fi

    # sort the SAM alignment file and convert it to BAM
    samtools sort \
      ${sample_id}_vs_${params.input.reference_prefix}.sam \
      -o ${sample_id}_vs_${params.input.reference_prefix}.bam \
      -O bam \
      -T temp

    # remove SAM file if it will not be published
    rm -f *.sam

    # index BAM alignment file
    samtools index ${sample_id}_vs_${params.input.reference_prefix}.bam
    samtools stats ${sample_id}_vs_${params.input.reference_prefix}.bam > ${sample_id}_vs_${params.input.reference_prefix}.bam.log

    # generate expression-level transcript abundance
    stringtie \
      -v \
      -p ${task.cpus} \
      -e \
      -o ${sample_id}_vs_${params.input.reference_prefix}.gtf \
      -G ${gtf_file} \
      -A ${sample_id}_vs_${params.input.reference_prefix}.ga \
      -l ${sample_id} ${sample_id}_vs_${params.input.reference_prefix}.bam

    # remove BAM file if it will not be published
    if [[ ${params.output.publish_bam} == false ]]; then
      rm -f *.bam
    fi

    # generate raw counts from hisat2/stringtie
    # Run the prepDE.py script provided by stringtie to get the raw counts.
    echo "${sample_id}\t./${sample_id}_vs_${params.input.reference_prefix}.gtf" > gtf_files
    prepDE.py -i gtf_files -g ${sample_id}_vs_${params.input.reference_prefix}.raw.pre

    # Reformat the raw file to be the same as the TPM/FKPM files.
    cat ${sample_id}_vs_${params.input.reference_prefix}.raw.pre | \
      grep -v gene_id | \
      perl -pi -e "s/,/\\t/g" > ${sample_id}_vs_${params.input.reference_prefix}.raw

    # generate the final FPKM and TPM files
    if [[ ${params.output.publish_fpkm} == true ]]; then
      awk -F"\t" '{if (NR!=1) {print \$1, \$8}}' OFS='\t' ${sample_id}_vs_${params.input.reference_prefix}.ga > ${sample_id}_vs_${params.input.reference_prefix}.fpkm
    fi

    if [[ ${params.output.publish_tpm} == true ]]; then
      awk -F"\t" '{if (NR!=1) {print \$1, \$9}}' OFS='\t' ${sample_id}_vs_${params.input.reference_prefix}.ga > ${sample_id}_vs_${params.input.reference_prefix}.tpm
    fi

  # or use kallisto
  elif [[ ${params.software.alignment} == 1 ]]; then
    # perform Kallisto alignment of fastq files
    if [ -e ${sample_id}_2.fastq ]; then
      kallisto quant \
        -i  ${params.input.reference_prefix}.transcripts.Kallisto.indexed \
        -o ${sample_id}_vs_${params.input.reference_prefix}.ga \
        ${sample_id}_1.fastq \
        ${sample_id}_2.fastq > ${sample_id}.kallisto.log 2>&1
    else
      kallisto quant \
        --single \
        -l 70 \
        -s .0000001 \
        -i ${params.input.reference_prefix}.transcripts.Kallisto.indexed \
        -o ${sample_id}_vs_${params.input.reference_prefix}.ga \
        ${sample_id}_1.fastq > ${sample_id}.kallisto.log 2>&1
    fi

    # generate TPM and raw count files
    if [[ ${params.output.publish_tpm} == true ]]; then
      awk -F"\t" '{if (NR!=1) {print \$1, \$5}}' OFS='\t' ${sample_id}_vs_${params.input.reference_prefix}.ga/abundance.tsv > ${sample_id}_vs_${params.input.reference_prefix}.tpm
    fi

    if [[ ${params.output.publish_raw} == true ]]; then
      awk -F"\t" '{if (NR!=1) {print \$1, \$4}}' OFS='\t' ${sample_id}_vs_${params.input.reference_prefix}.ga/abundance.tsv > ${sample_id}_vs_${params.input.reference_prefix}.raw
    fi

  # or use salmon
  elif [[ ${params.software.alignment} == 2 ]]; then
    # perform SALMON alignment of fastq files
    if [ -e ${sample_id}_2.fastq ]; then
      salmon quant \
        -i . \
        -l A \
        -1 ${sample_id}_1.fastq \
        -2 ${sample_id}_2.fastq \
        -p ${task.cpus} \
        -o ${sample_id}_vs_${params.input.reference_prefix}.ga \
        --minAssignedFrags 1 > ${sample_id}.salmon.log 2>&1
    else
      salmon quant \
        -i . \
        -l A \
        -r ${sample_id}_1.fastq \
        -p ${task.cpus} \
        -o ${sample_id}_vs_${params.input.reference_prefix}.ga \
        --minAssignedFrags 1 > ${sample_id}.salmon.log 2>&1
    fi

    # generate final TPM and raw count files
    if [[ ${params.output.publish_tpm} == true ]]; then
      awk -F"\t" '{if (NR!=1) {print \$1, \$4}}' OFS='\t' ${sample_id}_vs_${params.input.reference_prefix}.ga/quant.sf > ${sample_id}_vs_${params.input.reference_prefix}.tpm
    fi

    if [[ ${params.output.publish_raw} == true ]]; then
      awk -F"\t" '{if (NR!=1) {print \$1, \$5}}' OFS='\t' ${sample_id}_vs_${params.input.reference_prefix}.ga/quant.sf > ${sample_id}_vs_${params.input.reference_prefix}.raw
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
  # FPKM format is only generated if hisat2 is used
  if [[ ${params.output.publish_fpkm} == true && ${params.software.alignment} == 0 ]]; then
    create-gem.py \
      --sources ${workflow.launchDir}/${params.output.dir} \
      --prefix ${params.project.machine_name} \
      --type FPKM
  fi;

  if [[ ${params.output.publish_raw} == true ]]; then
    create-gem.py \
      --sources ${workflow.launchDir}/${params.output.dir} \
      --prefix ${params.project.machine_name} \
      --type raw
  fi

  if [[ ${params.output.publish_tpm} == true ]]; then
    create-gem.py \
      --sources ${workflow.launchDir}/${params.output.dir} \
      --prefix ${params.project.machine_name} \
      --type TPM
  fi
  """
}
