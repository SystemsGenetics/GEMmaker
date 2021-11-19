// Import generic module functions
include { saveFiles } from './functions'

params.options = [:]

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
    container "systemsgenetics/gemmaker:2.0.0"

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

    trimmomatic.sh \
        ${sample_id} \
        ${params.trimmomatic_MINLEN} \
        ${task.cpus} \
        ${fasta_adapter} \
        ${params.trimmomatic_LEADING} \
        ${params.trimmomatic_TRAILING} \
        ${params.trimmomatic_SLIDINGWINDOW} \
        "${fastq_files}"
    """
}
