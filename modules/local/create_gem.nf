// Import generic module functions
include { saveFiles } from './functions'

params.options = [:]

/**
 * Creates the GEM file from all the FPKM/TPM outputs
 */
process create_gem {
    publishDir "${params.outdir}/GEMs", mode: params.publish_dir_mode
    container "systemsgenetics/gemmaker:2.0.0"

    input:
    val(signal)
    path(input_files)

    output:
    path("*.GEM.*.txt"), emit: GEM_FILES

    script:
    """
    echo "#TRACE publish_fpkm=${params.publish_fpkm}"
    echo "#TRACE hisat2_enable=${params.hisat2_enable}"
    echo "#TRACE publish_tpm=${params.publish_tpm}"
    echo "#TRACE publish_raw=${params.publish_raw}"
    echo "#TRACE fpkm_lines=`cat ${params.outdir}/*.fpkm 2> /dev/null  | wc -l`"
    echo "#TRACE tpm_lines=`cat ${params.outdir}/*.tpm 2> /dev/null | wc -l`"
    echo "#TRACE raw_lines=`cat ${params.outdir}/*.raw 2> /dev/null | wc -l`"

    create_gem.sh \
        ${params.publish_fpkm} \
        ${params.hisat2_enable} \
        . \
        GEMmaker \
        ${params.publish_raw} \
        ${params.publish_tpm}
    """
}
