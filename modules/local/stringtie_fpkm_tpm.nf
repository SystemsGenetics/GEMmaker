/**
 * Generates the final FPKM / TPM / raw files from Hisat2
 */
process stringtie_fpkm_tpm {
    tag { sample_id }
    publishDir "${params.outdir}/Samples/${sample_id}", mode: params.publish_dir_mode
    // Note: 2.7X indices incompatible with AWS iGenomes.
    conda     (params.enable_conda ? "bioconda::stringtie=2.1.7" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/stringtie:2.1.7--h978d192_0"
    } else {
        container "quay.io/biocontainers/stringtie:2.1.7--h978d192_0"
    }

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
    echo "#TRACE publish_fpkm=${params.keep_fpkm}"
    echo "#TRACE publish_tpm=${params.keep_tpm}"
    echo "#TRACE publish_raw=${params.keep_counts}"
    echo "#TRACE ga_lines=`cat *.ga | wc -l`"
    echo "#TRACE gtf_lines=`cat *.gtf | wc -l`"

    if [[ ${params.keep_fpkm} == true ]]; then
      awk -F"\t" '{if (NR!=1) {print \$1, \$8}}' OFS='\t' ${sample_id}.Hisat2.ga > ${sample_id}.Hisat2.fpkm
    fi

    if [[ ${params.keep_tpm} == true ]]; then
      awk -F"\t" '{if (NR!=1) {print \$1, \$9}}' OFS='\t' ${sample_id}.Hisat2.ga > ${sample_id}.Hisat2.tpm
    fi

    if [[ ${params.keep_counts} == true ]]; then
      # Run the prepDE.py script provided by stringtie to get the raw counts.
      echo -e "${sample_id}\t./${sample_id}.Hisat2.gtf" > gtf_files
      prepDE.py -i gtf_files -g ${sample_id}.raw.pre

      # Reformat the raw file to be the same as the TPM/FKPM files.
      cat ${sample_id}.raw.pre | \
        grep -v gene_id | \
        sed "s/,/\\t/g" > ${sample_id}.Hisat2.raw
    fi
    """
}
