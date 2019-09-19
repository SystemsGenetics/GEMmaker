#!/bin/bash

sample_id="$1"

#This script determines if there is 1 or 2 files present and concatonates them together
if ls *_1.fastq >/dev/null 2>&1; then
  cat *_1.fastq >> "${sample_id}_1.fastq"
fi

if ls *_2.fastq >/dev/null 2>&1; then
  cat *_2.fastq >> "${sample_id}_2.fastq"
fi
