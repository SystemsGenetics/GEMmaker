#!/bin/bash

sample_id="$1"
params_software_trimmomatic_MINLEN="$2"
task_cpus="$3"
fasta_adapter="$4"
params_software_trimmomatic_LEADING="$5"
params_software_trimmomatic_TRAILING="$6"
params_software_trimmomatic_SLIDINGWINDOW="$7"
params_software_trimmomatic_quality=" "

echo $params_software_trimmomatic_LEADING
echo $params_software_trimmomatic_quality

# This script calculates average length of fastq files.
total=0
# This if statement checks if the data is single or paired data, and checks length accordingly
# This script returns 1 number, which can be used for the minlen in trimmomatic
if [ -e ${sample_id}_1.fastq ] && [ -e ${sample_id}_2.fastq ]; then
  for fastq in ${sample_id}_1.fastq ${sample_id}_2.fastq; do
    a=`awk 'NR%4 == 2 {lengths[length($0)]++} END {for (l in lengths) {print l, lengths[l]}}' $fastq \
    | sort \
    | awk '{ print $0, $1*$2}' \
    | awk -v var="$params_software_trimmomatic_MINLEN" '{ SUM += $3 } { SUM2 += $2 } END { printf("%.0f", SUM / SUM2 * var)} '`
  total=($a + $total)
  done
  total=( $total / 2 )
  minlen=$total
elif [ -e ${sample_id}_1.fastq ]; then
  minlen=`awk 'NR%4 == 2 {lengths[length($0)]++} END {for (l in lengths) {print l, lengths[l]}}' ${sample_id}_1.fastq \
    | sort \
    | awk '{ print $0, $1*$2}' \
    | awk -v var="$params_software_trimmomatic_MINLEN" '{ SUM += $3 } { SUM2 += $2 } END { printf("%.0f", SUM / SUM2 * var)} '`
fi
if [ -e ${sample_id}_1.fastq ] && [ -e ${sample_id}_2.fastq ]; then
  java -Xmx512m org.usadellab.trimmomatic.Trimmomatic \
    PE \
    -threads ${task_cpus} \
    ${params_software_trimmomatic_quality} \
    ${sample_id}_1.fastq \
    ${sample_id}_2.fastq \
    ${sample_id}_1p_trim.fastq \
    ${sample_id}_1u_trim.fastq \
    ${sample_id}_2p_trim.fastq \
    ${sample_id}_2u_trim.fastq \
    ILLUMINACLIP:${fasta_adapter}:2:40:15 \
    LEADING:${params_software_trimmomatic_LEADING} \
    TRAILING:${params_software_trimmomatic_TRAILING} \
    SLIDINGWINDOW:${params_software_trimmomatic_SLIDINGWINDOW} \
    MINLEN:"$minlen" > ${sample_id}.trim.log 2>&1
else
  # For ease of the next steps, rename the reverse file to the forward.
  # since these are non-paired it really shouldn't matter.
  if [ -e ${sample_id}_2.fastq ]; then
    mv ${sample_id}_2.fastq ${sample_id}_1.fastq
  fi
  # Now run trimmomatic
  java -Xmx512m org.usadellab.trimmomatic.Trimmomatic \
    SE \
    -threads ${task_cpus} \
    ${params_software_trimmomatic_quality} \
    ${sample_id}_1.fastq \
    ${sample_id}_1u_trim.fastq \
    ILLUMINACLIP:${fasta_adapter}:2:40:15 \
    LEADING:${params_software_trimmomatic_LEADING} \
    TRAILING:${params_software_trimmomatic_TRAILING} \
    SLIDINGWINDOW:${params_software_trimmomatic_SLIDINGWINDOW} \
    MINLEN:"$minlen" > ${sample_id}.trim.log 2>&1
fi
