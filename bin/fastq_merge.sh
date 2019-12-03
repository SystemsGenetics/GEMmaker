#!/bin/bash

# This script determines if there is 1 or 2 files present and concatonates
# them together

sample_id="$1"
publish="$2"


# First, concatenate all of the set 1 files
files1=`ls *_1.fastq | grep -v ${sample_id} | sort`
for file in $files1; do
   echo "Concatenate file: ${file} to ${sample_id}_1.fastq"
   cat $file >> "${sample_id}_1.fastq"
done
echo "Done with ${sample_id}_1.fastq"

# Next, concatenate all of the set 2 files
files2=`ls *_2.fastq | grep -v ${sample_id} | sort`
for file in $files2; do
  echo "Concatenate file: ${file} to ${sample_id}_2.fastq"
  cat $file >> "${sample_id}_2.fastq"
done
echo "Done with ${sample_id}_2.fastq"
