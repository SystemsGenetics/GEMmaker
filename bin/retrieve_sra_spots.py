#!/usr/bin/env python

"""
Get the total number of spots in a sample from
the sample's metadata file.
"""
import argparse
import json
import retrieve_sra_metadata as rsm



if __name__ == "__main__":

    # parse command-line arguments
    parser = argparse.ArgumentParser()
    parser.add_argument("meta_dir", help="Base directory for SRA metadata")
    parser.add_argument("sample_id", help="Sample ID")

    args = parser.parse_args()

    # load metadata file
    sample_dir = rsm.get_accession_directory(args.meta_dir, args.sample_id)
    metadata = json.load(open("%s/%s.ncbi.meta.json" % (sample_dir, args.sample_id), "r"))

    # parse total_spots attribute from sample metadata
    print(metadata["Pool"]["Member"]["@spots"])
