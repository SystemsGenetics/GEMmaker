#!/bin/bash

sample_id="$1"
kallisto_index="$2"
params_input_reference_name="$3"

if [ -e ${sample_id}_2.fastq ]; then
  kallisto quant \
    -i ${kallisto_index} \
    -o ${sample_id}_vs_${params_input_reference_name}.Kallisto.ga \
    ${sample_id}_1.fastq \
    ${sample_id}_2.fastq > ${sample_id}.kallisto.log 2>&1
else
  kallisto quant \
    --single \
    -l 70 \
    -s .0000001 \
    -i ${kallisto_index} \
    -o ${sample_id}_vs_${params_input_reference_name}.Kallisto.ga \
    ${sample_id}_1.fastq > ${sample_id}.kallisto.log 2>&1
fi
