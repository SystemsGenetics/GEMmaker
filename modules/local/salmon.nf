// Import generic module functions
include { saveFiles } from './functions'

params.options = [:]

/**
 * Performs SALMON alignemnt of fastq files
 */
process salmon {
    tag { sample_id }
    publishDir "${params.outdir}/Samples/${sample_id}", mode: params.publish_dir_mode, pattern: params.publish_pattern_salmon_ga
    container "systemsgenetics/gemmaker:2.0.0"

    input:
    tuple val(sample_id), path(fastq_files)
    path(salmon_index)

    output:
    tuple val(sample_id), path("*.ga", type: "dir"), emit: GA_FILES
    tuple val(sample_id), path("*.ga", type: "dir"), emit: LOGS
    tuple val(sample_id), val(params.DONE_SENTINEL), emit: DONE_SIGNAL

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
