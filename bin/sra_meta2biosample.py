#!/usr/bin/env python3

"""
A Python script for converting GEMmaker sample JSON meta data to a tab delimited file

.. module:: GEMmaker
    :platform: UNIX, Linux
    :synopsis: This script recursively reads in all the JSON files in a sample directory
        and parses it to find NCBI BioSample details. Results are saved into a tab
        delimited file.
"""


import json
import glob
from pathlib import Path
import urllib
import xmltodict
import os
import argparse
import pandas as pd
import pprint
import json
import sys
import re

def parse_meta(meta):
   accession = None
   sample = {}

   # Get the run accession
   if '@accession' in meta.keys():
       accession = meta['@accession'] 
   else:
       return sample

   # If this is an SRA sample accession ID then continue
   regexp = re.compile(r'^.RS\d+$')
   if not regexp.search(accession):
       return sample

   if 'IDENTIFIERS' in meta.keys():
       if 'EXTERNAL_ID' in meta['IDENTIFIERS'].keys():
           external_ids = meta['IDENTIFIERS']['EXTERNAL_ID']
           if isinstance(external_ids, dict):
               external_ids = [external_ids]
           for external_id in external_ids:
               if ('@namespace' in external_id) & (external_id['@namespace'] == 'BioSample'):
                   sample['ncbi_biosample_accession'] = external_id['#text']

   # If there is no sample then skip this one.
   if not ('ncbi_biosample_accession' in sample.keys()):
       pp.pprint(meta)
       return sample

   # Add in any sample attributes
   if 'SAMPLE_ATTRIBUTES' in meta.keys():
       attrs = meta['SAMPLE_ATTRIBUTES']['SAMPLE_ATTRIBUTE']
       if isinstance(attrs, dict):
           attrs = [attrs]
       for attr in attrs: 
          sample[attr['TAG']] = attr['VALUE']

   return sample

    

if __name__ == "__main__":
    pp = pprint.PrettyPrinter(indent=4)

    # Parse command-line arguments
    parser = argparse.ArgumentParser()

    parser.add_argument("--dir", help="A directory containing JSON files. Use this argument as often as needed", required=True, action='append')
    parser.add_argument("--out", help="The output file name", required=True)
    args = parser.parse_args()

    samples = []

    num_files = 0;
    for metadir in args.dir:
       for path in Path(metadir).rglob('*.json'):
           num_files = num_files + 1
           meta_file = open(path)
           meta = json.load(meta_file)
           sample = parse_meta(meta)
           if sample:
               samples.append(sample)

    samples = pd.DataFrame(samples)
    total_found = samples.shape[0]
    samples.drop_duplicates(['ncbi_biosample_accession'],inplace=True)
    num_samples = samples.shape[0]
    print("Parsed {} files. Found {} samples with {} unique.".format(num_files, total_found, num_samples), file=sys.stderr)
    samples.to_csv(args.out, sep="\t", index=False)
