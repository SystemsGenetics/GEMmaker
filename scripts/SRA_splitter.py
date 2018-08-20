#!/usr/bin/env python
##
##  This is a helper script for splitting a very large file 
##  of NCBI SRAs into smaller chunk files for use with GEMmaker.
## 
##  Usage:
##    python SRA_splitter.py [file] [chunk_size]
## 
## The [file] argument should be the name of the file with the 
## list of SRA IDs. It should be a two-column tab-delimited file
## where the fist column is the SRA experiment ID and the second
## is the SRA run ID.
##
## This script will then split the file into files near the
## size specified by [chunk_size].  Because GEMmaker merges all
## runs for an experiment into a single sample, this script will
## ensure that each split file does not split runs belonging to the
## same experiment into two different files.
##
##

import pandas as pd
import sys
import glob, os

script, SRA_file, chunk_size = sys.argv
chunk_size = int(chunk_size)

# Remove existing SRA_IDs.*.txt files
for f in glob.glob("SRA_IDs*.txt"):
  os.remove(f)

# Import the samples and decide the chunks size
print("Reading file " + SRA_file)
samples = pd.read_table(SRA_file, sep='\t', names=['EXPID','RUNID'])

# Open the first chunk
print("Splitting into files with approximatly %d samples per file" % chunk_size)
chunk_num = 1
filename = "SRA_IDs.%02d.txt" % (chunk_num)
chunkf = open(filename, 'w')
print("Writing " + filename)
prev_exp = ''
increment_chunk = False
for index, sample in samples.iterrows():

  # Write the sample to the chunk file.
  curr_exp = sample[0]
  curr_run = sample[1]
  chunkf.write(curr_run + "\n")
  
  # Check if we have added chunk_size items to the
  # file. if so, then indicate we need to increment to the next
  # chunk.
  if ((index + 1) % chunk_size == 0):
    increment_chunk = True 

  # If we should increment to the next chunk file make sure that
  # we are not in the middle of an experiment. If so, wait until
  # the experiment ends.
  if (increment_chunk and curr_exp != prev_exp):
    chunkf.close()
    chunk_num = chunk_num + 1
    filename = "SRA_IDs.%02d.txt" % (chunk_num)
    chunkf = open(filename, 'w');
    print("Writing " + filename)
    increment_chunk = False

  # Set this experiment as the prev_exp 
  prev_exp = curr_exp

chunkf.close()
