#!/bin/bash

sample_id="$1"
kallisto_index="$2"
bootstrap_samples="$3"
task_cpus="$4"
fastq_files="$5"

# convert the incoming FASTQ file list to an array
read -a fastq_files <<< $fastq_files

if [ ${#fastq_files[@]} == 2 ]; then
  kallisto quant \
    -i ${kallisto_index} \
    -b ${bootstrap_samples} \
    -o ${sample_id}.Kallisto.ga \
    -t ${task_cpus} \
    ${fastq_files[0]} \
    ${fastq_files[1]} > ${sample_id}.kallisto.log 2>&1
else
  kallisto quant \
    --single \
    -l 70 \
    -s .0000001 \
    -i ${kallisto_index} \
    -b ${bootstrap_samples} \
    -o ${sample_id}.Kallisto.ga \
    -t ${task_cpus} \
    ${fastq_files[0]} > ${sample_id}.kallisto.log 2>&1
fi
