#!/bin/bash

sra=$1
params_ref_prefix=$2

awk -F"\t" '{if (NR!=1) {print $1, $8}}' OFS='\t' ${sra}_vs_${params_ref_prefix}.ga > ${sra}_vs_${params_ref_prefix}.fpkm
