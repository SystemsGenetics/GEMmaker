nextflow.enable.dsl=2

/**
 * Performs STAR alignment of fastq files to a genome reference
 */
process star {
    tag { sample_id }
    publishDir "${params.outdir}/Samples/${sample_id}", mode: params.publish_dir_mode, pattern: "*.log"

    conda (params.enable_conda ? "bioconda::star=2.7.9a" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
      container "https://depot.galaxyproject.org/singularity/star:2.7.9a--h9ee0642_0"
    } else {
      container "quay.io/biocontainers/star:2.7.9a--h9ee0642_0"
    }

    input:
    tuple val(sample_id), path(fastq_files)
    path(star_index)

    output:
    tuple val(sample_id), path("*.sam"), emit: SAM_FILES
    tuple val(sample_id), path("*Log.final.out"), emit: LOGS
    tuple val(sample_id), val(params.DONE_SENTINEL), emit: DONE_SIGNAL

    script:
    """
    echo "#TRACE sample_id=${sample_id}"
    echo "#TRACE n_cpus=${task.cpus}"
    echo "#TRACE trimmed_fastq_lines=`cat *.fastq | wc -l`"
    echo "#TRACE index_bytes=`stat -Lc '%s' ${star_index} | awk '{sum += \$1} END {print sum}'`"

    # convert the incoming FASTQ file list to an array
    fastq_files=(${fastq_files})

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
      # First align the paired.
      STAR \
        --runThreadN ${task.cpus} \
        --genomeDir ${star_index} \
        --outFileNamePrefix ${sample_id}.trimmed_paired \
        --readFilesIn \${fq_1p} \${fq_2p} > ${sample_id}.trimmed_paired.star.log 2>&1

      # Now aligned the non-paired
      STAR \
        --runThreadN ${task.cpus} \
        --genomeDir ${star_index} \
        --outFileNamePrefix ${sample_id}.trimmed_single \
        --readFilesIn \${fq_1u},\${fq_2u} > ${sample_id}.trimmed_single.star.log 2>&1

    else
        STAR \
            --runThreadN ${task.cpus} \
            --genomeDir ${star_index} \
            --outFileNamePrefix ${sample_id}.trimmed \
            --readFilesIn \${fq_1u} > ${sample_id}.trimmed.star.log 2>&1
    fi
    """
}
