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
    
    # FPKM format is only generated if hisat2/star is used
    if [[ ${params.publish_fpkm} == true ]]; then
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
