// Import generic module functions
include { saveFiles } from './functions'

params.options = [:]

/**
 * Generates the final FPKM / TPM / raw files from Hisat2
 */
process hisat2_fpkm_tpm {
    tag { sample_id }
    publishDir "${params.outdir}/Samples/${sample_id}", mode: params.publish_dir_mode
    container "systemsgenetics/gemmaker:2.0.0"

    input:
    tuple val(sample_id), path(stringtie_files)

    output:
    path("*.Hisat2.fpkm"), optional: true, emit: FPKM_FILES
    path("*.Hisat2.tpm"), optional: true, emit: TPM_FILES
    path("*.Hisat2.raw"), optional: true, emit: RAW_FILES
    tuple val(sample_id), val(params.DONE_SENTINEL), emit: DONE_SIGNAL

    script:
    """
    echo "#TRACE sample_id=${sample_id}"
    echo "#TRACE publish_fpkm=${params.hisat2_keep_fpkm}"
    echo "#TRACE publish_tpm=${params.hisat2_keep_tpm}"
    echo "#TRACE publish_raw=${params.hisat2_keep_counts}"
    echo "#TRACE ga_lines=`cat *.ga | wc -l`"
    echo "#TRACE gtf_lines=`cat *.gtf | wc -l`"

    hisat2_fpkm_tpm.sh \
        ${params.hisat2_keep_fpkm} \
        ${sample_id} \
        ${params.hisat2_keep_tpm} \
        ${params.hisat2_keep_counts}
    """
}
