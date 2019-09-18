#!/bin/bash

params_output_publish_tpm="$1"
sample_id="$2"
params_input_reference_name="$3"
params_output_publish_raw="$4"
  
if [[ ${params_output_publish_tpm} == true ]]; then
  awk -F"\t" '{if (NR!=1) {print $1, $4}}' OFS='\t' ${sample_id}_vs_${params_input_reference_name}.Salmon.ga/quant.sf > ${sample_id}_vs_${params_input_reference_name}.Salmon.tpm
fi

if [[ ${params_output_publish_raw} == true ]]; then
  awk -F"\t" '{if (NR!=1) {print $1, $5}}' OFS='\t' ${sample_id}_vs_${params_input_reference_name}.Salmon.ga/quant.sf > ${sample_id}_vs_${params_input_reference_name}.Salmon.raw
fi
