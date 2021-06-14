#!/bin/bash

sample_id="$1"
base_name="$2"
task_cpus="$3"
fastq_files="$4"

# convert the incoming FASTQ file list to an array
read -a fastq_files <<< $fastq_files

# we don't know the order the files will come so we have
# to find the paired and non paired files.
fq_1p=""
fq_2p=""
fq_1u=""
fq_2u=""
for f in "${fastq_files[@]}"; do
    echo $f
    if [[ $f =~ _1p_trim.fastq ]]; then
        fq_1p=$f
    elif [[ $f =~ _2p_trim.fastq ]]; then
        fq_2p=$f
    elif [[ $f =~ _1u_trim.fastq ]]; then
        fq_1u=$f
    elif [[ $f =~ _2u_trim.fastq ]]; then
        fq_2u=$f
    fi
done;

if [ ${#fastq_files[@]} == 4 ]; then
  hisat2 \
    -x ${base_name} \
    --no-spliced-alignment \
    -q \
    -1 ${fq_1p} \
    -2 ${fq_2p} \
    -U ${fq_1u},${fq_2u} \
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
    -U ${fq_1u} \
    -S ${sample_id}.sam \
    -t \
    -p ${task_cpus} \
    --un ${sample_id}_un.fastq \
    --dta-cufflinks \
    --new-summary \
    --summary-file ${sample_id}.sam.log
fi
