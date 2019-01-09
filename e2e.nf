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

General Information:
--------------------
  Profile(s):         ${workflow.profile}
  Container Engine:   ${workflow.containerEngine}


Input Parameters:
-----------------
  Remote fastq list path:     ${params.input.remote_list_path}
  Local sample glob:          ${params.input.local_samples_path}
  Reference genome path:      ${params.input.reference_path}
  Reference genome prefix:    ${params.input.reference_prefix}


Output Parameters:
------------------
  Output directory:           ${params.output.dir}
  Publish SRA:                ${params.output.publish_sra}
  Publish downloaded FASTQ:   ${params.output.publish_downloaded_fastq}
  Publish trimmed FASTQ:      ${params.output.publish_trimmed_fastq}
  Publish BAM:                ${params.output.publish_bam}
  Publish FPKM:               ${params.output.publish_fpkm}
  Publish TPM:                ${params.output.publish_tpm}


Execution Parameters:
---------------------
  Queue size:                 ${params.execution.queue_size}
  Number of threads:          ${params.execution.threads}
  Maximum retries:            ${params.execution.max_retries}
  Error strategy:             ${params.execution.error_strategy}


Software Parameters:
--------------------
  Trimmomatic clip path:      ${params.software.trimmomatic.clip_path}
  Trimmomatic minimum ratio:  ${params.software.trimmomatic.MINLEN}
"""



/**
 * Create value channels that can be reused
 */
HISAT2_INDEXES = Channel.fromPath("${params.input.reference_path}/${params.input.reference_prefix}*.ht2*").collect()
SALMON_INDEXES = Channel.fromPath("${params.input.reference_path}/${params.input.reference_prefix}.transcripts.Salmon.indexed/*").collect()
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
  SRR_FILE = Channel.value(params.input.remote_list_path)
}



/**
 * Retrieves metadata for all of the remote samples
 * and maps SRA runs to SRA experiments.
 */
process retrieve_sample_metadata {
  publishDir params.output.dir, mode: params.output.publish_mode, pattern: "*.GEMmaker.meta.*", saveAs: { "${it.tokenize(".")[0]}/${it}" }
  label "python3"

  input:
    val srr_file from SRR_FILE

  output:
    stdout REMOTE_SAMPLES_LIST
    file "*.GEMmaker.meta.*"

  script:
  """
  retrieve_SRA_metadata.py $srr_file
  """
}



/**
 * Create channels for remote samples and local samples
 */
REMOTE_SAMPLES_LIST
  .splitCsv()
  .groupTuple(by: 1)
  .map { [it[1], it[0].toString().replaceAll(/[\[\]\'\,]/,''), 'remote'] }
  .set{REMOTE_SAMPLES}

LOCAL_SAMPLE_FILES
  .map{ [it[0], it[1], 'local' ] }
  .set{LOCAL_SAMPLES}

// ALL_SAMPLES = REMOTE_SAMPLES.mix(LOCAL_SAMPLES)



if ( params.output.publish_fpkm == false && params.output.publish_tpm == false ) {
  error "At least one of FPKM or TPM should be enabled."
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
  label "retry"
  publishDir params.output.sample_dir, mode: params.output.publish_mode

  input:
    set val(sample_id), val(run_ids), val(type) from REMOTE_SAMPLES
    file indexes from HISAT2_INDEXES
    file gtf_file from GTF_FILE

  output:
    val(sample_id) into COMPLETED_SAMPLES
    file("*.sra") optional true into SRA_FILES
    file("*.fastq") optional true into FASTQ_FILES
    file("*fastqc.*") into FASTQC_FILES
    file("*.log") into LOG_FILES
    file("*.sam") optional true into SAM_FILES
    file("*.bam") optional true into BAM_FILES
    file("*.bam.bai") optional true into BAI_FILES
    file("*.ga") into GA_FILES
    file("*.gtf") into GTF_FILES
    file("*.raw") into RAW_FILES
    file("*.fpkm") optional true into FPKM_FILES
    file("*.tpm") optional true into TPM_FILES

  script:
  """
  # download SRA files from NCBI
  SRR_IDS=\$(echo ${run_ids} | perl -p -e 's/[\\[,\\]]//g')

  # use ascp
  if [[ ${params.software.sra_download} == 0 ]]; then
    for id in \$SRR_IDS; do
      prefix=`echo \$id | perl -p -e 's/^([SDE]RR)\\d+\$/\$1/'`
      sixchars=`echo \$id | perl -p -e 's/^([SDE]RR\\d{1,3}).*\$/\$1/'`
      url="/sra/sra-instant/reads/ByRun/sra/\$prefix/\$sixchars/\$id/\$id.sra"
      ascp -i /home/gemdocker/.aspera/connect/etc/asperaweb_id_dsa.openssh -k 1 -T -l 1000m anonftp@ftp.ncbi.nlm.nih.gov:\$url .
    done

  # or use the SRA toolkit
  elif [[ ${params.software.sra_download} == 1 ]]; then
    for id in \$SRR_IDS; do
      prefetch --output-directory . \$id
    done
  fi

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

  # TODO: run the above code only for remote samples

  # perform fastqc on raw FASTQ files
  MERGED_FASTQ_FILES=\$(ls ${sample_id}_?.fastq)

  fastqc \$MERGED_FASTQ_FILES

  # TODO: salmon, kallisto

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
      -threads ${params.execution.threads} \
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
      -threads ${params.execution.threads} \
      ${params.software.trimmomatic.quality} \
      ${sample_id}_1.fastq \
      ${sample_id}_1u_trim.fastq \
      ILLUMINACLIP:${params.software.trimmomatic.clip_path}:2:40:15 \
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
      -p ${params.execution.threads} \
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
      -p ${params.execution.threads} \
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
  if [[ ${params.output.publish_sam} == false ]]; then
    rm -f *.sam
  fi

  # index BAM alignment file
  samtools index ${sample_id}_vs_${params.input.reference_prefix}.bam
  samtools stats ${sample_id}_vs_${params.input.reference_prefix}.bam > ${sample_id}_vs_${params.input.reference_prefix}.bam.log

  # generate expression-level transcript abundance
  stringtie \
    -v \
    -p ${params.execution.threads} \
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

  script:
  """
  multiqc --ignore ${params.output.dir}/GEM --ignore ${params.output.dir}/reports ${params.output.dir}
  """
}



/**
 * Creates the GEM file from all the FPKM/TPM outputs
 */
process create_gem {
  label "python3"
  publishDir "${params.output.dir}/GEM", mode: params.output.publish_mode

  input:
    val signal from COMPLETED_SAMPLES_FOR_GEM.collect()

  output:
    file "*.GEM.*.txt" into GEM_FILES

  script:
  """
  # If the alignment tool is hisat then we need to generate both
  # TPM and FPKM
  if [ ${params.software.alignment.which_alignment} == 0 ]; then
    create_GEM.py --sources ${params.output.dir} --prefix ${params.project.machine_name} --type FPKM
  fi;
  create_GEM.py --sources ${params.output.dir} --prefix ${params.project.machine_name} --type raw
  create_GEM.py --sources ${params.output.dir} --prefix ${params.project.machine_name} --type TPM
  """
}
