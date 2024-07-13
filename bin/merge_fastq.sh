#!/bin/bash

sample_id=$1

files_1=(`find . -maxdepth 1 -name "*_1.fastq" | grep -v $sample_id | sort`)
files_2=(`find . -maxdepth 1 -name "*_2.fastq" | grep -v $sample_id | sort`)

# Check the number of files. If there is more than one then
# combine them. If not, then just create a link with the
# sample_id as the name.
if [ ${#files_1[@]} -gt 1 ] ; then

    for i in "${!files_1[@]}"; do
        cat ${files_1[$i]} >> ${sample_id}_1.fastq
    done;

elif [ ${#files_1[@]} -eq 1 ] ; then

    ln ${files_1[0]} ${sample_id}_1.fastq

fi

# If there are paired files.
if [ ${#files_2[@]} -gt 1 ] ; then

    for i in "${!files_2[@]}"; do
        cat ${files_2[$i]} >> ${sample_id}_2.fastq
    done;

elif [ ${#files_2[@]} -eq 1 ] ; then

    ln ${files_2[0]} ${sample_id}_2.fastq
fi

