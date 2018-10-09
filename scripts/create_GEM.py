#!/usr/bin/env python

"""A Python script for creating gene expression matrix (GEM) files


.. module:: GEMmaker
  :platform: Unix, Windows
  :synopsis: The Gene expression matrix (GEM) file can be created after a successful run of 
    GEMmaker.  It can be used to create GEM files containing either FPKM or TPM
    values.  It can also be used to combine the results from multiple GEMmaker 
    directories into a single GEM file.

.. moduleauthor:: John Hadish & Stephen Ficklin

"""
import argparse
import pandas as pd
import glob
import os
import re

# Specify the arguments that are allowed by this script
parser = argparse.ArgumentParser(description='Process some integers.')
parser.add_argument('--sources', dest='path', action='store', required=True, nargs='*',
                   help='One or more GEMmaker directory paths where results are stored. FPKM or TPM files are found in Sample_* or [SDR]RX directories within the GEMmaker directory.')
parser.add_argument('--prefix', dest='prefix', action='store', required=True,
                   help='A prefix to be added to the GEM.txt file when writing the final matrix.')
parser.add_argument('--type', dest='type', action='store', required=True, choices=['TPM', 'FPKM'],
                   help='The type of count values to include in the  GEM (either "TPM" or "FPKM".')

# Read in the input arguments
args = parser.parse_args()

# Set the name of the output GEM file.
ematrix_name = args.prefix + ".GEM." + args.type + '.txt';

# Iterate through the GEMmaker directories and find the FPKM and TPM files
result_files = []
for source_dir in args.path:
  # Remove any trailing slash from the source dir.
  source_dir = re.sub(r'/$', '', source_dir)

  # Check for NCBI SRA files
  print("Finding %s files in %s" % (args.type,  source_dir))
  for filename in glob.iglob(source_dir + '/[SED]RX*/*' + str.lower(args.type), recursive=True):
    result_files.append(filename)

  # Check for local non SRA files
  for filename in glob.iglob(source_dir + '/Sample_*/*' + str.lower(args.type), recursive=True):
    result_files.append(filename)

print("Found %d sample files" % (len(result_files)))

# Initialize the expression matrix by reading in the
ematrix = pd.DataFrame({'gene' : []})
for result in result_files:

  # Get the sample name from the file name.
  file_basename = os.path.basename(result)
  sample_name = file_basename.split('_vs_')[0]

  # Remove the Sample_ from local filenames.
  sample_name = re.sub(r'^Sample_', r'', sample_name)
  print ("Adding results for sample: "  + sample_name)
  df = pd.read_csv(result, header = None, sep = '\t', names = ["gene", sample_name])

  # Stringtie was creating duplicate entries for some genes because of 
  # how it handles genes that have non-overlapping splice variants. So, we need
  # to remove the duplicates but keep the first occurance.
  df.drop_duplicates(['gene'], keep = 'first', inplace = True)

  # Now merge in this sample to the expression matrix.
  ematrix = ematrix.merge(df.iloc[:,[0,1]], on='gene', how='outer')

# Set the gene names as the data frame indexes.
ematrix = ematrix.set_index('gene', drop=True)

# Write out our expression matrix
print("Writing " + ematrix_name)
ematrix.to_csv(ematrix_name, sep = '\t', na_rep="NA", index_label=False)

