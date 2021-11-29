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
    publishDir "${params.outdir}/Samples/${sample_id}", mode: params.publish_dir_mode, pattern: params.publish_pattern_trimmomatic
    container "systemsgenetics/gemmaker:2.1.0"

    input:
    tuple val(sample_id), path(fastq_files)
    path(fasta_adapter)

    output:
    tuple val(sample_id), path("*_trim.fastq"), emit: FASTQ_FILES
    tuple val(sample_id), path("*.trim.log"), emit: LOGS
    tuple val(sample_id), val(params.DONE_SENTINEL), emit: DONE_SIGNAL

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

    # convert the incoming FASTQ file list to an array
    read -a fastq_files <<< ${fastq_files}

    # This script calculates average length of fastq files.
    total=0
    # This if statement checks if the data is single or paired data, and checks length accordingly
    # This script returns 1 number, which can be used for the minlen in trimmomatic
    if [ \${#fastq_files[@]} == 2 ]; then
      for fastq in "\${fastq_files[@]}"; do
        cat="cat \$fastq"
        if [[ \$fastq =~ .gz\$ ]]; then
          cat="zcat \$fastq"
        fi
        a=`\$cat | awk 'NR%4 == 2 {lengths[length(\$0)]++} END {for (l in lengths) {print l, lengths[l]}}' \
        | sort \
        | awk '{ print \$0, \$1*\$2}' \
        | awk -v var="${params.trimmomatic_MINLEN}" '{ SUM += \$3 } { SUM2 += \$2 } END { printf("%.0f", SUM / SUM2 * var)} '`
      total=(\$a + \$total)
      done
      total=( \$total / 2 )
      minlen=\$total
    else
      cat="cat \${fastq_files[0]}"
      if [[ \${fastq_files[0]} =~ .gz\$ ]]; then
        cat="zcat \${fastq_files[0]}"
      fi
      minlen=`\$cat | awk 'NR%4 == 2 {lengths[length(\$0)]++} END {for (l in lengths) {print l, lengths[l]}}'  \
        | sort \
        | awk '{ print \$0, \$1*\$2}' \
        | awk -v var="${params.trimmomatic_MINLEN}" '{ SUM += \$3 } { SUM2 += \$2 } END { printf("%.0f", SUM / SUM2 * var)} '`
    fi
    if [ \${#fastq_files[@]} == 2 ]; then
      trimmomatic \
        PE \
        -threads ${task.cpus} \
        \${fastq_files[0]} \
        \${fastq_files[1]} \
        ${sample_id}_1p_trim.fastq \
        ${sample_id}_1u_trim.fastq \
        ${sample_id}_2p_trim.fastq \
        ${sample_id}_2u_trim.fastq \
        ILLUMINACLIP:${fasta_adapter}:2:40:15 \
        LEADING:${params.trimmomatic_LEADING} \
        TRAILING:${params.trimmomatic_TRAILING} \
        SLIDINGWINDOW:${params.trimmomatic_SLIDINGWINDOW} \
        MINLEN:"\$minlen" > ${sample_id}.trim.log 2>&1
    else
      trimmomatic \
        SE \
        -threads ${task.cpus} \
        \${fastq_files[0]} \
        ${sample_id}_1u_trim.fastq \
        ILLUMINACLIP:${fasta_adapter}:2:40:15 \
        LEADING:${params.trimmomatic_LEADING} \
        TRAILING:${params.trimmomatic_TRAILING} \
        SLIDINGWINDOW:${params.trimmomatic_SLIDINGWINDOW} \
        MINLEN:"\$minlen" > ${sample_id}.trim.log 2>&1
    fi
    """
}
