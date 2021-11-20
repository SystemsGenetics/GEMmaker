/**
 * Process to generate the multiqc report once everything is completed
 */
process multiqc {
    cache false
    publishDir "${params.outdir}/reports", mode: params.publish_dir_mode
    
    conda (params.enable_conda ? 'bioconda::multiqc=1.11' : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/multiqc:1.11--pyhdfd78af_0"
    } else {
        container "quay.io/biocontainers/multiqc:1.11--pyhdfd78af_0"
    }

    input:
    val(signal)
    path(config_file)
    path(custom_logo)
    path(input_files)

    output:
    path("multiqc_data"), emit: MULTIQC_DATA
    path("multiqc_report.html"), emit: MULTIQC_REPORT

    script:
    """
    multiqc --config ${config_file} ./
    """
}
