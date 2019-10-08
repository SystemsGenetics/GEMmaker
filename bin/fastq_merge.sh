#!/bin/bash

# This script determines if there is 1 or 2 files present and concatonates
# them together

sample_id="$1"
publish="$2"

# Cleanup any previously existing sample FASTQ files
if [ -e ${sample_id}_1.fastq ] && [ ! -e ${sample_id}_1.done ]; then
  rm ${sample_id}_1.fastq
fi
if [ -e ${sample_id}_2.fastq ] && [ ! -e ${sample_id}_2.done ]; then
  rm ${sample_id}_2.fastq
fi

# First, concatenate all of the set 1 files
if [ ! -e ${sample_id}_1.done ]; then
  files1=`ls *_1.fastq | grep -v ${sample_id} | sort`
  for file in $files1; do
     cat $file >> "${sample_id}_1.fastq"
  done
  touch ${sample_id}_1.done
fi

# Next, concatenate all of the set 2 files 
if [ ! -e ${sample_id}_2.done ]; then
   files2=`ls *_2.fastq | grep -v ${sample_id} | sort`
   for file in $files2; do
      # Copy the contents of the run FASTQ into the sample FASTQ
      cat $file >> "${sample_id}_2.fastq"
   done
   touch ${sample_id}_2.done
fi

# If we're all done, then clean away the original files
if [ -e ${sample_id}_1.done ] && [ -e ${sample_id}_2.done ]; then
   files=`ls *.fastq | grep -v ${sample_id} | sort`
   for file in $files; do 
     # Get the linked path for the original run FASTQ file and
     # cleanup the run files.
     lpath=`stat -c %N  $file | awk -F"'" '{print $4}'`
     echo "clean_work_files.sh $lpath"
     if [ $publish == "null" ]; then
       clean_work_files.sh $lpath
     fi
  done
fi
