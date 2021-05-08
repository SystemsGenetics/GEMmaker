#!/bin/bash

sample_id="$1"
task_cpus="$2"
salmon_index="$3"
fastq_files="$4"

# convert the incoming FASTQ file list to an array
read -a fastq_files <<< $fastq_files

if [ ${#fastq_files[@]} == 2 ]; then
  salmon quant \
    -i ${salmon_index} \
    -l A \
    -1 ${fastq_files[0]} \
    -2 ${fastq_files[1]} \
    -p ${task_cpus} \
    -o ${sample_id}.Salmon.ga \
    --minAssignedFrags 1 > ${sample_id}.salmon.log 2>&1
else
  salmon quant \
    -i ${salmon_index} \
    -l A \
    -r ${fastq_files[0]} \
    -p ${task_cpus} \
    -o ${sample_id}.Salmon.ga \
    --minAssignedFrags 1 > ${sample_id}.salmon.log 2>&1
fi

# Copy these files for MultiQC reporting
cp ${sample_id}.Salmon.ga/aux_info/meta_info.json ${sample_id}-meta_info.json
cp ${sample_id}.Salmon.ga/libParams/flenDist.txt ${sample_id}-flenDist.txt
