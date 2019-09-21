#!/bin/bash

params_output_publish_fpkm="$1"
sample_id="$2"
params_input_reference_name="$3"
params_output_publish_tpm="$4"
params_output_publish_raw="$5"
  
if [[ ${params_output_publish_fpkm} == true ]]; then
  awk -F"\t" '{if (NR!=1) {print $1, $8}}' OFS='\t' ${sample_id}_vs_${params_input_reference_name}.Hisat2.ga > ${sample_id}_vs_${params_input_reference_name}.Hisat2.fpkm
fi

if [[ ${params_output_publish_tpm} == true ]]; then
  awk -F"\t" '{if (NR!=1) {print $1, $9}}' OFS='\t' ${sample_id}_vs_${params_input_reference_name}.Hisat2.ga > ${sample_id}_vs_${params_input_reference_name}.Hisat2.tpm
fi

if [[ ${params_output_publish_raw} == true ]]; then
  # Run the prepDE.py script provided by stringtie to get the raw counts.
  echo -e "${sample_id}\t./${sample_id}_vs_${params_input_reference_name}.Hisat2.gtf" > gtf_files
  prepDE.py -i gtf_files -g ${sample_id}_vs_${params_input_reference_name}.raw.pre

  # Reformat the raw file to be the same as the TPM/FKPM files.
  cat ${sample_id}_vs_${params_input_reference_name}.raw.pre | \
    grep -v gene_id | \
    perl -pi -e "s/,/\\t/g" > ${sample_id}_vs_${params_input_reference_name}.Hisat2.raw
  fi
