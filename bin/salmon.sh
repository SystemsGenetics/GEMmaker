#!/bin/bash

sample_id="$1"
task_cpus="$2"

if [ -e ${sample_id}_1.fastq ] && [ -e ${sample_id}_2.fastq ]; then
  salmon quant \
    -i . \
    -l A \
    -1 ${sample_id}_1.fastq \
    -2 ${sample_id}_2.fastq \
    -p ${task_cpus} \
    -o ${sample_id}.Salmon.ga \
    --minAssignedFrags 1 > ${sample_id}.salmon.log 2>&1
elif [ -e ${sample_id}_1.fastq ]; then
  salmon quant \
    -i . \
    -l A \
    -r ${sample_id}_1.fastq \
    -p ${task_cpus} \
    -o ${sample_id}.Salmon.ga \
    --minAssignedFrags 1 > ${sample_id}.salmon.log 2>&1
else
  salmon quant \
    -i . \
    -l A \
    -r ${sample_id}_2.fastq \
    -p ${task_cpus} \
    -o ${sample_id}.Salmon.ga \
    --minAssignedFrags 1 > ${sample_id}.salmon.log 2>&1
fi
