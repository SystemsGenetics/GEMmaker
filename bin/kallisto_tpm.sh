#!/bin/bash

sample_id="$1"
params_output_publish_tpm="$2"
params_output_publish_raw="$3"
params_input_reference_name="$4"

if [[ ${params_output_publish_tpm} == true ]]; then
  awk -F"\t" '{if (NR!=1) {print $1, $5}}' OFS='\t' ${sample_id}_vs_${params_input_reference_name}.Kallisto.ga/abundance.tsv > ${sample_id}_vs_${params_input_reference_name}.Kallisto.tpm
fi

if [[ ${params_output_publish_raw} == true ]]; then
  awk -F"\t" '{if (NR!=1) {print $1, $4}}' OFS='\t' ${sample_id}_vs_${params_input_reference_name}.Kallisto.ga/abundance.tsv > ${sample_id}_vs_${params_input_reference_name}.Kallisto.raw
fi
