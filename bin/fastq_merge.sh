#!/bin/bash

# This script determines if there is 1 or 2 files present and concatonates
# them together

sample_id="$1"
publish="$2"

# Cleanup any previously existing sample FASTQ files
if [ -e ${sample_id}_1.fastq ]; then
  rm ${sample_id}_1.fastq
fi
if [ -e ${sample_id}_2.fastq ]; then
  rm ${sample_id}_2.fastq
fi

files1=`ls *_1.fastq | grep -v ${sample_id} | sort`
for file in $files1; do
   # Copy the contents of the run FASTQ into the sample FASTQ
   cat $file >> "${sample_id}_1.fastq"
   
   # Get the linked path for the original run FASTQ file and
   # cleanup the run files.
   lpath=`stat -c %N  $file | awk -F"'" '{print $4}'`
   echo "clean_work_files.sh $lpath"
   if [ $publish == "null" ]; then
     clean_work_files.sh $lpath
   fi
done

files2=`ls *_2.fastq | grep -v ${sample_id} | sort`
for file in $files2; do
   # Copy the contents of the run FASTQ into the sample FASTQ
   cat $file >> "${sample_id}_2.fastq"
   
   # Get the linked path for the original run FASTQ file and
   # cleanup the run files.
   lpath=`stat -c %N  $file | awk -F"'" '{print $4}'`
   if [ $publish == "null" ]; then
      clean_work_files.sh $lpath
   fi
done
