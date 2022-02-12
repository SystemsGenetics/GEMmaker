/**
 * Sorts the SAM alignment file (or BAM file in the case of STAR paired end reads) and coverts it to binary BAM
 */
process samtools_sort {
    tag { sample_id }
    publishDir "${params.outdir}/Samples/${sample_id}", mode: params.publish_dir_mode, pattern: params.publish_pattern_samtools_sort

    conda (params.enable_conda ? "bioconda::samtools=1.14" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/samtools:1.14--hb421002_0"
    } else {
        container "quay.io/biocontainers/samtools:1.14--hb421002_0"
    }

    input:
    tuple val(sample_id), path(sam_file)

    output:
    tuple val(sample_id), path("*sorted.bam"), emit: BAM_FILES
    tuple val(sample_id), val(params.DONE_SENTINEL), emit: DONE_SIGNAL

    script:
    """
    echo "#TRACE sample_id=${sample_id}"
    # echo "#TRACE sam_lines=`cat *.sam | wc -l`"

    samtools sort \
        -o ${sample_id}.sorted.bam \
        -O bam \
        -T temp \
        ${sam_file}
    """
}
