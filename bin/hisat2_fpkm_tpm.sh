#!/bin/bash

params_output_keep_fpkm="$1"
sample_id="$2"
params_output_keep_tpm="$3"
params_output_keep_counts="$4"

if [[ ${params_output_keep_fpkm} == true ]]; then
  awk -F"\t" '{if (NR!=1) {print $1, $8}}' OFS='\t' ${sample_id}.Hisat2.ga > ${sample_id}.Hisat2.fpkm
fi

if [[ ${params_output_keep_tpm} == true ]]; then
  awk -F"\t" '{if (NR!=1) {print $1, $9}}' OFS='\t' ${sample_id}.Hisat2.ga > ${sample_id}.Hisat2.tpm
fi

if [[ ${params_output_keep_counts} == true ]]; then
  # Run the prepDE.py script provided by stringtie to get the raw counts.
  echo -e "${sample_id}\t./${sample_id}.Hisat2.gtf" > gtf_files
  prepDE.py -i gtf_files -g ${sample_id}.raw.pre

  # Reformat the raw file to be the same as the TPM/FKPM files.
  cat ${sample_id}.raw.pre | \
    grep -v gene_id | \
    perl -pi -e "s/,/\\t/g" > ${sample_id}.Hisat2.raw
fi
