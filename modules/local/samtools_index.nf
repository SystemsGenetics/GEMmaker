/**
 * Indexes the BAM alignment file
 */
process samtools_index {
    tag { sample_id }
    publishDir "${params.outdir}/Samples/${sample_id}", mode: params.publish_dir_mode, pattern: params.publish_pattern_samtools_index
    
    conda (params.enable_conda ? "bioconda::samtools=1.14" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/samtools:1.14--hb421002_0"
    } else {
        container "quay.io/biocontainers/samtools:1.14--hb421002_0"
    }

    input:
    tuple val(sample_id), path(bam_file)

    output:
    tuple val(sample_id), path(bam_file), emit: BAM_FILES
    tuple val(sample_id), path("*.bam.bai"), emit: BAI_FILES
    tuple val(sample_id), path("*.bam.log"), emit: LOGS

    script:
    """
    echo "#TRACE sample_id=${sample_id}"
    echo "#TRACE bam_bytes=`stat -Lc '%s' *.bam`"

    samtools index ${bam_file}
    samtools stats ${bam_file} > ${sample_id}.bam.log
    """
}
