// Import generic module functions
include { saveFiles } from './functions'

params.options = [:]

/**
 * Generates the final TPM and raw count files for Kallisto
 */
process kallisto_tpm {
    tag { sample_id }
    publishDir "${params.outdir}/Samples/${sample_id}", mode: params.publish_dir_mode
    container "systemsgenetics/gemmaker:2.0.0"

    input:
    tuple val(sample_id), path(ga_file)

    output:
    path("*.tpm"), optional: true, emit: TPM_FILES
    path("*.raw"), optional: true, emit: RAW_FILES
    tuple val(sample_id), val(params.DONE_SENTINEL), emit: DONE_SIGNAL

    script:
    """
    echo "#TRACE sample_id=${sample_id}"

    kallisto_tpm.sh \
        ${sample_id} \
        ${params.kallisto_keep_tpm} \
        ${params.kallisto_keep_counts}
    """
}
