/**
 * Creates the GEM file from all the FPKM/TPM outputs
 */
process create_gem {
    publishDir "${params.outdir}/GEMs", mode: params.publish_dir_mode
    container "systemsgenetics/gemmaker:2.1.0"

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

    # FPKM format is only generated if hisat2 is used
    if [[ ${params.publish_fpkm} == true && ${params.hisat2_enable} == true ]]; then
      create-gem.py \
        --sources . \
        --prefix GEMmaker \
        --type FPKM
    fi;

    if [[ ${params.publish_raw} == true ]]; then
      create-gem.py \
        --sources . \
        --prefix GEMmaker \
        --type raw
    fi

    if [[ ${params.publish_tpm} == true ]]; then
      create-gem.py \
        --sources . \
        --prefix GEMmaker \
        --type TPM
    fi
    """
}
