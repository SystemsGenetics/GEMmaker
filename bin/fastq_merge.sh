#!/bin/bash

# This script determines if there is 1 or 2 files present and concatonates
# them together

sample_id="$1"
publish="$2"

files1=`ls *_1.fastq | sort`

for file in $files1; do
   cat $file >> "${sample_id}_1.fastq"
   lpath=`stat -c %N  $file | awk -F"'" '{print $4}'`
   echo "clean_work_files.sh $lpath"
   if [ $publish == "null" ]; then
     clean_work_files.sh $lpath
   fi
done

files2=`ls *_2.fastq | sort`
for file in $files2; do
   cat $file >> "${sample_id}_2.fastq"
   lpath=`stat -c %N  $file | awk -F"'" '{print $4}'`
   if [ $publish == "null" ]; then
      clean_work_files.sh $lpath
   fi
done
