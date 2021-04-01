#!/bin/bash
# Run pipeline script with the example data

nextflow run main.nf -ansi-log false -profile travis
