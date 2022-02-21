#!/usr/bin/env python3

"""
A Python script for converting NCBI JSON meta data to a tab delimited file
.. module:: GEMmaker
    :platform: UNIX, Linux
    :synopsis: This script recursively reads in all the JSON files in a single directory
        and parses it to find NCBI BioSample details. Results are saved into a tab
        delimited file
"""

import json
import glob
from pathlib import Path
import os
import argparse
import pandas as pd
import pprint
import json


def parse_meta(meta):
   sample = {}
   if 'IDENTIFIERS' in meta.keys():
       if 'EXTERNAL_ID' in meta['IDENTIFIERS'].keys():
           external_id = meta['IDENTIFIERS']['EXTERNAL_ID']
           if ('@namespace' in external_id) & (external_id['@namespace'] == 'BioSample'):
               sample['ncbi_biosample_accession'] = external_id['#text']

   # If there is no sample then skip this one.
   if not ('ncbi_biosample_accession' in sample.keys()):
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

    parser.add_argument("--dir", help="A directory containing JSON files", required=True)
    parser.add_argument("--out", help="The output file name", required=True)
    args = parser.parse_args()

    samples = []

    for path in Path(args.dir).rglob('*.json'):
       meta_file = open(path)
       meta = json.load(meta_file)
       try:
           sample = parse_meta(meta)
           samples.append(sample)
           #pp.pprint(sample)
       except:
         pp.pprint(meta)

    samples = pd.DataFrame(samples)
    print(samples)
    samples.to_csv(args.out, sep="\t", index=False)
