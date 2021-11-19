// Import generic module functions
include { saveFiles } from './functions'

params.options = [:]

/**
 * Indexes the BAM alignment file
 */
process samtools_index {
    tag { sample_id }
    publishDir "${params.outdir}/Samples/${sample_id}", mode: params.publish_dir_mode, pattern: params.publish_pattern_samtools_index
    container "systemsgenetics/gemmaker:2.0.0"

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
