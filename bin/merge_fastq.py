#!/usr/bin/env python3

"""A Python script for merging multiple FASTQ files


.. module:: GEMmaker
  :platform: Unix, Windows
  :synopsis: This merges multiple FASTQ files into a single file. It will
  keep paired files separate. Output files will be named {out_prefix}_1.fastq
  and {out_prefix}_2.fastq, qhere {out_prefix} is the value of the
  --out_prefix argument.

.. moduleauthor:: Stephen Ficklin

"""
from Bio import SeqIO
import sys
import argparse
import re

parser = argparse.ArgumentParser(description='Process some integers.')
parser.add_argument("--fastq_files", dest='fastq_files', type=str, required=True,
                    help="List of input FASTQ files", nargs="+")
parser.add_argument("--out_prefix", dest='out_prefix', type=str, required=True,
                    help="The prefix for the ouptput FASTQ files.")
args = parser.parse_args()

# Keep an array of file handles based on pairs.
handles = {
 '_1.fastq': None,
 '_2.fastq': None
}

# Only open file handles for the paired file suffix we find.
for fastq_filename in args.fastq_files:
    if (re.match(r"^.*(_1.fastq)$", fastq_filename) and handles['_1.fastq'] is None):
        handles['_1.fastq'] = open(args.out_prefix + '_1.fastq', "w")
    if (re.match(r"^.*(_2.fastq)$", fastq_filename) and handles['_2.fastq'] is None):
        handles['_2.fastq'] = open(args.out_prefix + '_2.fastq', "w")

# Now open each file and merge the records into the appropirate
# output file.
for fastq_filename in args.fastq_files:
    fh = None
    if (re.match(r"^.*(_1.fastq)$", fastq_filename)):
        fh = handles['_1.fastq']
    if (re.match(r"^.*(_2.fastq)$", fastq_filename)):
        fh = handles['_2.fastq']
    with open(fastq_filename) as fastq_file:
        for record in SeqIO.parse(fastq_file, "fastq"):
            SeqIO.write(record, fh, "fastq")

# Close the output file handles.
if (handles['_1.fastq'] is not None):
    handles['_1.fastq'].close()
if (handles['_2.fastq'] is not None):
    handles['_2.fastq'].close()
