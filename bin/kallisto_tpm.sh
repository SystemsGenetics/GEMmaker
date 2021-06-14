#!/bin/bash

sample_id="$1"
params_output_keep_tpm="$2"
params_output_keep_counts="$3"

if [[ ${params_output_keep_tpm} == true ]]; then
  awk -F"\t" '{if (NR!=1) {print $1, $5}}' OFS='\t' ${sample_id}.Kallisto.ga/abundance.tsv > ${sample_id}.Kallisto.tpm
fi

if [[ ${params_output_keep_counts} == true ]]; then
  awk -F"\t" '{if (NR!=1) {print $1, $4}}' OFS='\t' ${sample_id}.Kallisto.ga/abundance.tsv > ${sample_id}.Kallisto.raw
fi
