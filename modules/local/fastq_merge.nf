// Import generic module functions
include { saveFiles } from './functions'

params.options = [:]

/**
 * This process merges the fastq files based on their sample_id number.
 */
process fastq_merge {
    tag { sample_id }
    container "systemsgenetics/gemmaker:2.0.0"

    input:
    tuple val(sample_id), path(fastq_files)

    output:
    tuple val(sample_id), path("${sample_id}_?.fastq"), emit: FASTQ_FILES
    tuple val(sample_id), val(params.DONE_SENTINEL), emit: DONE_SIGNAL

    script:
    """
    echo "#TRACE sample_id=${sample_id}"
    echo "#TRACE fastq_lines=`cat *.fastq | wc -l`"

    fastq_merge.sh ${sample_id}
    """
}
