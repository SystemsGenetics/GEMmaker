#!/usr/bin/env python
##
## This is a helper script for splitting a very large file
## of NCBI SRAs into smaller chunk files for use with GEMmaker.
##
## Usage:
##   python SRA_splitter.py [file] [chunk-size]
##
## The [file] argument should be the name of the file with the
## list of SRA IDs. It should be a two-column tab-delimited file
## where the fist column is the SRA experiment ID and the second
## is the SRA run ID.
##
## This script will then split the file into files with at most
## [chunk-size] lines per file. Because GEMmaker merges all
## runs for an experiment into a single sample, this script ensures
## that the runs of an experiment are not split between two files.
##
import glob
import os
import pandas as pd
import sys

# parse command-line arguments
if len(sys.argv) != 3:
  print("usage: python SRA_splitter.py [file] [chunk-size]")

SRA_file = sys.argv[1]
chunk_size = int(sys.argv[2])

# remove existing SRA_IDs.*.txt files
for f in glob.glob("SRA_IDs.*.txt"):
  os.remove(f)

# load the SRA file as a dataframe
print("Reading file %s" % (SRA_file))

samples = pd.read_table(SRA_file, names=["exp_id", "run_id"])

# open the first chunk file
chunk_num = 1
chunk_filename = "SRA_IDs.%02d.txt" % (chunk_num)
chunk_file = open(chunk_filename, "w")
curr_size = 0

# iterate through unique experiment IDs
experiments = list(set(samples["exp_id"]))
experiments.sort()

for exp in experiments:
  # fetch the run IDs for the given experiment ID
  runs = samples.loc[samples["exp_id"] == exp, "run_id"]

  # throw error if there are too many run IDs for chunk size
  if len(runs) > chunk_size:
    print("error: experiment %s has %d runs which cannot be satisfied by chunk size of %d" % (exp, len(runs), chunk_size))
    sys.exit(-1)

  # determine whether the current chunk file has enough
  # lines remaining for the given run IDs
  if curr_size + len(runs) > chunk_size:
    # close the current chunk file
    print("Wrote %4d runs to %s" % (curr_size, chunk_filename))
    chunk_file.close()

    # open the next chunk file
    chunk_num += 1
    chunk_filename = "SRA_IDs.%02d.txt" % (chunk_num)
    chunk_file = open(chunk_filename, "w")
    curr_size = 0

  # write the run IDs to the chunk file
  chunk_file.write("".join(["%s\t%s\n" % (exp, run) for run in runs]))
  curr_size += len(runs)

# close the last chunk file
print("Wrote %4d runs to %s" % (curr_size, chunk_filename))
chunk_file.close()
