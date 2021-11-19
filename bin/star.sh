#!/bin/bash

sample_id="$1"
genome_dir="$2"
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
  # First align the paired.
  STAR \
    --runThreadN ${task_cpus} \
    --genomeDir ${genome_dir} \
    --readFilesIn ${fq_1p} ${fq_2p}
  # Now aligned the non-paired
  STAR \
    --runThreadN ${task_cpus} \
    --genomeDir ${genome_dir} \
    --readFilesIn ${fq_1u},${fq_2u}

else
    STAR \
        --runThreadN ${task_cpus} \
        --genomeDir ${genome_dir} \
        --readFilesIn ${fq_1u}
fi
