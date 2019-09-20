#!/bin/bash


params_output_publish_fpkm="$1"
params_software_alignment="$2"
workflow_launchDir="$3"
params_output_dir="$4"
params_project_machine_name="$5"
params_output_publish_raw="$6"
params_output_publish_tpm="$7"

# FPKM format is only generated if hisat2 is used
if [[ ${params_output_publish_fpkm} == true && ${params_software_alignment} == 0 ]]; then
  create-gem.py \
    --sources ${workflow_launchDir}/${params_output_dir} \
    --prefix ${params_project_machine_name} \
    --type FPKM
fi;

if [[ ${params_output_publish_raw} == true ]]; then
  create-gem.py \
    --sources ${workflow_launchDir}/${params_output_dir} \
    --prefix ${params_project_machine_name} \
    --type raw
fi

if [[ ${params_output_publish_tpm} == true ]]; then
  create-gem.py \
    --sources ${workflow_launchDir}/${params_output_dir} \
    --prefix ${params_project_machine_name} \
    --type TPM
fi
