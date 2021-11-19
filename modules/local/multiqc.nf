// Import generic module functions
include { saveFiles } from './functions'

params.options = [:]

/**
 * Process to generate the multiqc report once everything is completed
 */
process multiqc {
    cache false
    publishDir "${params.outdir}/reports", mode: params.publish_dir_mode
    container "systemsgenetics/gemmaker:2.0.0"

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
