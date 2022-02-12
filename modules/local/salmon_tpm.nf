/**
 * Generates the final TPM file for Salmon
 */
process salmon_tpm {
    tag { sample_id }
    publishDir "${params.outdir}/Samples/${sample_id}", mode: params.publish_dir_mode
    container "systemsgenetics/gemmaker:2.1.0"

    input:
    tuple val(sample_id), path(ga_file)

    output:
    path("*.Salmon.tpm"), optional: true, emit: TPM_FILES
    path("*.Salmon.raw"), optional: true, emit: RAW_FILES
    tuple val(sample_id), val(params.DONE_SENTINEL), emit: DONE_SIGNAL

    script:
    """
    echo "#TRACE sample_id=${sample_id}"

    if [[ ${params.salmon_keep_tpm} == true ]]; then
      awk -F"\t" '{if (NR!=1) {print \$1, \$4}}' OFS='\t' ${sample_id}.Salmon.ga/quant.sf > ${sample_id}.Salmon.tpm
    fi

    if [[ ${params.salmon_keep_counts} == true ]]; then
      awk -F"\t" '{if (NR!=1) {print \$1, \$5}}' OFS='\t' ${sample_id}.Salmon.ga/quant.sf > ${sample_id}.Salmon.raw
    fi
    """
}
