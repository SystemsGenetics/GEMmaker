// Import generic module functions
include { saveFiles } from './functions'

params.options = [:]

/**
 * Performs hisat2 alignment of fastq files to a genome reference
 */
process star {
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
    tuple val(sample_id), val(params.DONE_SENTINEL), emit: DONE_SIGNAL

    script:
    """
    echo "#TRACE sample_id=${sample_id}"
    echo "#TRACE n_cpus=${task.cpus}"
    echo "#TRACE trimmed_fastq_lines=`cat *.fastq | wc -l`"
    echo "#TRACE index_bytes=`stat -Lc '%s' ${indexes} | awk '{sum += \$1} END {print sum}'`"

    star.sh \
        ${sample_id} \
        ${star_index} \
        ${task.cpus} \
        "${fastq_files}" \
    """
}
