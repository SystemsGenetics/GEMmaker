#!/bin/bash

sample_id="$1"
base_name="$2"
task_cpus="$3"


if [ -e ${sample_id}_2p_trim.fastq ]; then
  hisat2 \
    -x ${base_name} \
    --no-spliced-alignment \
    -q \
    -1 ${sample_id}_1p_trim.fastq \
    -2 ${sample_id}_2p_trim.fastq \
    -U ${sample_id}_1u_trim.fastq,${sample_id}_2u_trim.fastq \
    -S ${sample_id}.sam \
    -t \
    -p ${task_cpus} \
    --un ${sample_id}_un.fastq \
    --dta-cufflinks \
    --new-summary \
    --summary-file ${sample_id}.sam.log
else
  hisat2 \
    -x ${base_name} \
    --no-spliced-alignment \
    -q \
    -U ${sample_id}_1u_trim.fastq \
    -S ${sample_id}.sam \
    -t \
    -p ${task_cpus} \
    --un ${sample_id}_un.fastq \
    --dta-cufflinks \
    --new-summary \
    --summary-file ${sample_id}.sam.log
fi
