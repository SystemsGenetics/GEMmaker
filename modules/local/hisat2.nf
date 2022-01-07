/**
 * Performs hisat2 alignment of fastq files to a genome reference
 */
process hisat2 {
    tag { sample_id }
    publishDir "${params.outdir}/Samples/${sample_id}", mode: params.publish_dir_mode, pattern: "*.log"

    conda (params.enable_conda ? "bioconda::hisat2=2.2.0 bioconda::samtools=1.10" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/mulled-v2-a97e90b3b802d1da3d6958e0867610c718cb5eb1:2880dd9d8ad0a7b221d4eacda9a818e92983128d-0"
    } else {
        container "quay.io/biocontainers/mulled-v2-a97e90b3b802d1da3d6958e0867610c718cb5eb1:2880dd9d8ad0a7b221d4eacda9a818e92983128d-0"
    }

    input:
    tuple val(sample_id), path(fastq_files)
    path(indexes)

    output:
    tuple val(sample_id), path("*.sam"), emit: SAM_FILES
    tuple val(sample_id), path("*.sam.log"), emit: LOGS
    tuple val(sample_id), val(params.DONE_SENTINEL), emit: DONE_SIGNAL

    script:
    """
    echo "#TRACE sample_id=${sample_id}"
    echo "#TRACE n_cpus=${task.cpus}"
    echo "#TRACE trimmed_fastq_lines=`cat *.fastq | wc -l`"
    echo "#TRACE index_bytes=`stat -Lc '%s' ${indexes} | awk '{sum += \$1} END {print sum}'`"

    # convert the incoming FASTQ file list to an array
    read -a fastq_files <<< ${fastq_files}

    # we don't know the order the files will come so we have
    # to find the paired and non paired files.
    fq_1p=""
    fq_2p=""
    fq_1u=""
    fq_2u=""
    for f in "\${fastq_files[@]}"; do
        echo \$f
        if [[ \$f =~ _1p_trim.fastq ]]; then
            fq_1p=\$f
        elif [[ \$f =~ _2p_trim.fastq ]]; then
            fq_2p=\$f
        elif [[ \$f =~ _1u_trim.fastq ]]; then
            fq_1u=\$f
        elif [[ \$f =~ _2u_trim.fastq ]]; then
            fq_2u=\$f
        fi
    done;

    if [ \${#fastq_files[@]} == 4 ]; then
      hisat2 \
        -x ${params.hisat2_base_name} \
        --no-spliced-alignment \
        -q \
        -1 \${fq_1p} \
        -2 \${fq_2p} \
        -U \${fq_1u},\${fq_2u} \
        -S ${sample_id}.sam \
        -t \
        -p ${task.cpus} \
        --un ${sample_id}_un.fastq \
        --dta-cufflinks \
        --new-summary \
        --summary-file ${sample_id}.sam.log
    else
      hisat2 \
        -x ${params.hisat2_base_name} \
        --no-spliced-alignment \
        -q \
        -U \${fq_1u} \
        -S ${sample_id}.sam \
        -t \
        -p ${task.cpus} \
        --un ${sample_id}_un.fastq \
        --dta-cufflinks \
        --new-summary \
        --summary-file ${sample_id}.sam.log
    fi
    """
}
