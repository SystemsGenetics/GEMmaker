/**
 * Generates expression-level transcript abundance
 */
process stringtie {
    tag { sample_id }
    label "process_medium"
    publishDir "${params.outdir}/Samples/${sample_id}", mode: params.publish_dir_mode, pattern: params.publish_pattern_stringtie_ga_gtf

    // Note: 2.7X indices incompatible with AWS iGenomes.
    conda     (params.enable_conda ? "bioconda::stringtie=2.1.7" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/stringtie:2.1.7--h978d192_0"
    } else {
        container "quay.io/biocontainers/stringtie:2.1.7--h978d192_0"
    }

    input:
    tuple val(sample_id), path(bam_file)
    path(gtf_file)

    output:
    tuple val(sample_id), path("*.stringtie.*"), emit: GA_GTF_FILES
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
        -o ${sample_id}.stringtie.gtf \
        -G ${gtf_file} \
        -A ${sample_id}.stringtie.ga \
        -l ${sample_id} \
        ${bam_file}
    """
}
