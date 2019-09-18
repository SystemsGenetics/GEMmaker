#!/bin/bash

sample_id="$1"
params_input_hisat2_index_prefix="$2"
params_input_reference_name="$3"
task_cpus="$4"


if [ -e ${sample_id}_2p_trim.fastq ]; then
  hisat2 \
    -x ${params_input_hisat2_index_prefix} \
    --no-spliced-alignment \
    -q \
    -1 ${sample_id}_1p_trim.fastq \
    -2 ${sample_id}_2p_trim.fastq \
    -U ${sample_id}_1u_trim.fastq,${sample_id}_2u_trim.fastq \
    -S ${sample_id}_vs_${params_input_reference_name}.sam \
    -t \
    -p ${task_cpus} \
    --un ${sample_id}_un.fastq \
    --dta-cufflinks \
    --new-summary \
    --summary-file ${sample_id}_vs_${params_input_reference_name}.sam.log
else
  hisat2 \
    -x ${params_input_hisat2_index_prefix} \
    --no-spliced-alignment \
    -q \
    -U ${sample_id}_1u_trim.fastq \
    -S ${sample_id}_vs_${params_input_reference_name}.sam \
    -t \
    -p ${task_cpus} \
    --un ${sample_id}_un.fastq \
    --dta-cufflinks \
    --new-summary \
    --summary-file ${sample_id}_vs_${params_input_reference_name}.sam.log
fi
