#!/bin/bash

sra=$1
params_ref_prefix=$2

awk -F"\t" '{if ($3 == "transcript") print $0}' ${sra}_vs_${params_ref_prefix}.gtf | perl -p -e 's/^.*?transcript_id "(.*?)";.*FPKM "(.    *?)";.*$/$1\t$2/' > ${sra}_vs_${params_ref_prefix}.fpkm

