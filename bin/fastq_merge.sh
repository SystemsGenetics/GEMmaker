#!/bin/bash

sample_id="$1"

#This script determines if there is 1 or 2 files present and concatonates them together
if ls *_1.fastq >/dev/null 2>&1; then
  number_1_files=`ls *_1.fastq |wc -w` #This determines the number of files.
  if [ $number_1_files -gt 1 ]; then # If there are more than 1 files ending in _1, cat them into 1 file.

  # Look at the following beautiful bash code!!! I probably should have gone to python, but just look at how beautiful it is!
  # This determines the largest file, renames it to the proper output file, then concatonates everything onto the
  # end of the output. As you can tell, I am very happy with my bash skills.
    mv `ls -S *_1.fastq|head -n1` "${sample_id}_1.fastq"
    find . -maxdepth 1 -iname '*_1.fastq' -not -name "${sample_id}_1.fastq" -exec cat {} +>>"${sample_id}_1.fastq"
    echo "Option 1.1: Concatonate more than 1 '*_1' files" #This can be removed, but it is not harming the program. It is used for testing.
  else
    mv *_1.fastq "${sample_id}_1.fastq" # If there is only 1 file, rename it to the sample ID
    echo "Option 1.2: Only one '*_1' file, just rename" #This can be removed, but it is not harming the program. It is used for testing.
  fi
fi
#Same thing here except for _2
if ls *_2.fastq >/dev/null 2>&1; then
  number_2_files=`ls *_2.fastq |wc -w`
  if [ $number_2_files -gt 1 ]; then
    mv `ls -S *_2.fastq|head -n1` "${sample_id}_2.fastq"
    find . -maxdepth 1 -iname '*_2.fastq' -not -name "${sample_id}_2.fastq" -exec cat {} +>>"${sample_id}_2.fastq"
    echo "Option 2.1: Concatonate more than 1 '*_2' files" #This can be removed, but it is not harming the program. It is used for testing.
  else
    mv *_2.fastq "${sample_id}_2.fastq"
    echo "Option 2.2: Only one '*_2' file, just rename" #This can be removed, but it is not harming the program. It is used for testing.
  fi
fi
