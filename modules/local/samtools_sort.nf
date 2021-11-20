/**
 * Sorts the SAM alignment file and coverts it to binary BAM
 */
process samtools_sort {
    tag { sample_id }
    publishDir "${params.outdir}/Samples/${sample_id}", mode: params.publish_dir_mode, pattern: params.publish_pattern_samtools_sort
    container "systemsgenetics/gemmaker:2.0.0"

    input:
    tuple val(sample_id), path(sam_file)

    output:
    tuple val(sample_id), path("*.bam"), emit: BAM_FILES
    tuple val(sample_id), val(params.DONE_SENTINEL), emit: DONE_SIGNAL

    script:
    """
    echo "#TRACE sample_id=${sample_id}"
    echo "#TRACE sam_lines=`cat *.sam | wc -l`"

    samtools sort \
        -o ${sample_id}.bam \
        -O bam \
        -T temp \
        ${sam_file}
    """
}
