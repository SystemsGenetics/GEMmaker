#!/bin/bash

set -e

# make sure input directory doesn't already exist
if [[ -d input ]]; then
	echo "error: input directory already exists"
	exit -1
fi

# copy example data to input directory
cp -r examples/CORG input

# run each pipeline script with the example data
nextflow run main.nf -ansi-log false -profile travis
nextflow run main-e2e.nf -ansi-log false -profile travis
