#!/bin/bash

params_output_keep_tpm="$1"
params_output_keep_counts="$2"
sample_id="$3"

if [[ ${params_output_keep_tpm} == true ]]; then
  awk -F"\t" '{if (NR!=1) {print $1, $4}}' OFS='\t' ${sample_id}.Salmon.ga/quant.sf > ${sample_id}.Salmon.tpm
fi

if [[ ${params_output_keep_counts} == true ]]; then
  awk -F"\t" '{if (NR!=1) {print $1, $5}}' OFS='\t' ${sample_id}.Salmon.ga/quant.sf > ${sample_id}.Salmon.raw
fi
