// Import generic module functions
include { saveFiles } from './functions'

params.options = [:]

/**
 * Performs KALLISTO alignemnt of fastq files
 */
process kallisto {
    tag { sample_id }
    publishDir "${params.outdir}/Samples/${sample_id}", mode: params.publish_dir_mode, pattern: params.publish_pattern_kallisto_ga
    container "systemsgenetics/gemmaker:2.0.0"

    input:
    tuple val(sample_id), path(fastq_files)
    path(kallisto_index)

    output:
    tuple val(sample_id), path("*.ga", type: "dir"), emit: GA_FILES
    tuple val(sample_id), path("*.kallisto.log"), emit: LOGS
    tuple val(sample_id), val(params.DONE_SENTINEL), emit: DONE_SIGNAL

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
