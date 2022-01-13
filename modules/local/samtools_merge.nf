/**
 * Merge Multiple output files from STAR into 1. This is because we are running Trimmomatic
 */
process samtools_merge {
    tag { sample_id }
    publishDir "${params.outdir}/Samples/${sample_id}", mode: params.publish_dir_mode, pattern: params.publish_pattern_samtools_merge

    conda (params.enable_conda ? "bioconda::samtools=1.14" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/samtools:1.14--hb421002_0"
    } else {
        container "quay.io/biocontainers/samtools:1.14--hb421002_0"
    }

    input:
    tuple val(sample_id), path(sam_files)

    output:
    tuple val(sample_id), path("*.bam"), emit: BAM_FILES
    tuple val(sample_id), val(params.DONE_SENTINEL), emit: DONE_SIGNAL

    script:
    """
    echo "#TRACE sample_id=${sample_id}"
    echo "#TRACE sam_lines=`cat *.sam | wc -l`"

    samtools merge \
        -o ${sample_id}.bam \
        ${sam_files}
    """
}
