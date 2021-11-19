// Import generic module functions
include { saveFiles } from './functions'

params.options = [:]


/**
 * Generates expression-level transcript abundance
 */
process stringtie {
    tag { sample_id }
    label "process_medium"
    publishDir "${params.outdir}/Samples/${sample_id}", mode: params.publish_dir_mode, pattern: params.publish_pattern_stringtie_ga_gtf
    container "systemsgenetics/gemmaker:2.0.0"

    input:
    tuple val(sample_id), path(bam_file)
    path(gtf_file)

    output:
    tuple val(sample_id), path("*.Hisat2.*"), emit: GA_GTF_FILES
    tuple val(sample_id), val(params.DONE_SENTINEL), emit: DONE_SIGNAL

    script:
    """
    echo "#TRACE sample_id=${sample_id}"
    echo "#TRACE bam_bytes=`stat -Lc '%s' *.bam`"
    echo "#TRACE gtf_lines=`cat *.gtf | wc -l`"

    stringtie \
        -v \
        -p ${task.cpus} \
        -e \
        -o ${sample_id}.Hisat2.gtf \
        -G ${gtf_file} \
        -A ${sample_id}.Hisat2.ga \
        -l ${sample_id} \
        ${bam_file}
    """
}
