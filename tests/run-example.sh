#!/bin/bash

set -e

# run each pipeline script with the example data
nextflow run main.nf -ansi-log false -profile travis
nextflow run main-e2e.nf -ansi-log false -profile travis
