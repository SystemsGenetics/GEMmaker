/**
 * Performs fastqc on raw fastq files
 */
process fastqc {
    tag { sample_id }
    publishDir "${params.outdir}/Samples/${sample_id}", mode: params.publish_dir_mode, pattern: "*_fastqc.*"
    container "systemsgenetics/gemmaker:2.0.0"

    input:
    tuple val(sample_id), path(fastq_files)

    output:
    tuple path("*_fastqc.html"), path("*_fastqc.zip"), emit: REPORTS
    tuple val(sample_id), val(params.DONE_SENTINEL), emit: DONE_SIGNAL

    script:
    """
    echo "#TRACE sample_id=${sample_id}"
    echo "#TRACE fastq_lines=`cat *.fastq | wc -l`"

    fastqc ${fastq_files}
    """
}
