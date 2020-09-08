#!/bin/bash

# This script determines if there is 1 or 2 files present and concatonates
# them together. It also renames a single sample if the first pair is missing

sample_id="$1"


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


# If there is a FASTQ sample with _2 suffix but no _1  then rename
# to _1 so that count software will work
if [ -e ${sample_id}_2 ] && [ ! -e ${sample_id}_2 ]; then
  mv ${sample_id}_2 ${sample_id}_1
fi
