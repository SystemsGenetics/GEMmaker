/**
 * Creates a report of any SRA run IDs that failed and why they failed.
 */
process failed_run_report {
    publishDir "${params.outdir}/reports", mode: params.publish_dir_mode
    container "systemsgenetics/gemmaker:2.1.0"

    input:
    path(failed_runs)
    path(failed_run_template)

    output:
    path("failed_SRA_run_report.html"), emit: FAILED_RUN_REPORT

    script:
    """
    failed_runs_report.py --template ${failed_run_template}
    """
}
