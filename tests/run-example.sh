#!/bin/bash

set -e

# create input directory
mkdir input

# copy example data to input directory
cp -r examples/reference input
cp examples/LocalRunExample/**/*.fastq input
cp examples/RemoteRunExample/SRA_IDs.txt input

# run each pipeline script with the example data
nextflow run main.nf -ansi-log false -profile travis
nextflow run main-e2e.nf -ansi-log false -profile travis
