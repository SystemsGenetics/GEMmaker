#!/bin/bash


params_output_publish_fpkm="$1"
params_input_hisat2_enable="$2"
params_output_dir="$3"
params_project_machine_name="$4"
params_output_publish_raw="$5"
params_output_publish_tpm="$6"

# FPKM format is only generated if hisat2 is used
if [[ ${params_output_publish_fpkm} == true && ${params_input_hisat2_enable} == true ]]; then
  create-gem.py \
    --sources ${params_output_dir} \
    --prefix ${params_project_machine_name} \
    --type FPKM
fi;

if [[ ${params_output_publish_raw} == true ]]; then
  echo "create-gem.py --sources ${params_output_dir}  --prefix ${params_project_machine_name} --type raw"
  create-gem.py \
    --sources ${params_output_dir} \
    --prefix ${params_project_machine_name} \
    --type raw
fi

if [[ ${params_output_publish_tpm} == true ]]; then
  create-gem.py \
    --sources ${params_output_dir} \
    --prefix ${params_project_machine_name} \
    --type TPM
fi
