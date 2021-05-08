#!/bin/bash

sample_id="$1"
kallisto_index="$2"
task_cpus="$3"
fastq_files="$4"

# convert the incoming FASTQ file list to an array
read -a fastq_files <<< $fastq_files

if [ ${#fastq_files[@]} == 2 ]; then
  kallisto quant \
    -i ${kallisto_index} \
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
    -o ${sample_id}.Kallisto.ga \
    -t ${task_cpus} \
    ${fastq_files[0]} > ${sample_id}.kallisto.log 2>&1
fi
