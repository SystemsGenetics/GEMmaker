#!/usr/bin/env python3

"""
A Python script for creating an annotation matrix for NCBI experiments in a GEM.

.. module:: GEMmaker
    :platform: UNIX, Linux
    :synopsis: This script recieves a single input argument: a file containing
       a list of run IDs (numbers with SRR, ERR or DRR prefixes) from NCBI's
       sequence read archive (SRA). It generates a variety of JSON meta files
       and outputs a tab-delimited file that maps run IDs to experiment IDs.
"""
import argparse
import codecs
import json
import pandas as pd
import pprint
import re
import sys
import urllib
import xmltodict
import pprint
import pathlib
import os





def map_annotations(experiment, sample, meta_dir, cv_map):
    """
    Maps NCBI SRA annotations to controlled vocabularies.

    This function creates both a JSON and tab delimited metadata file that
    maps metadata from NCBI to known controlled vocabulary terms. The purpose
    of this is to help ensure uniformity in metadata details between
    runs of GEMmaker and between different sample sets.

    :param experiment: A dictionary containing the experiment metadata.
    :param sample: A dictionary containing the sample metadata.
    :param run: A dictionary containing the run metatdata.

    :return: a dictionary of exerpiment annotations
    """

    # For case-insenitivity, convert all SRA tags to lower-case
    tags = [x.lower() for x in cv_map.index.values]

    # Now that we have all the metadata loaded create the non-nested
    # annotation dictionary
    annots = {}

    # Build info about the biological sample.
    # Term: biological sample, sep:00195
    # Term: data:2091, Accession
    # Term: schema:title
    # Term: schema:name
    annots["sep:00195"] = {}
    annots["sep:00195"]["data:2091"] = sample.get("@accession", "")
    annots["sep:00195"]["schema:title"] = sample.get("TITLE", "")
    annots["sep:00195"]["schema:name"] = sample.get("@alias", "")

    # Add the organism and it's child terms.
    # Term: organism, obi:0100026
    # Term: rdfs:label
    # Term: Scientific Name, NCIT:C43459
    # Term: NCBI Taxonomy ID, data:1179
    annots["sep:00195"]["obi:0100026"] = {}
    annots["sep:00195"]["obi:0100026"]["rdfs:label"] = sample["SAMPLE_NAME"].get("SCIENTIFIC_NAME", "")
    annots["sep:00195"]["obi:0100026"]["NCIT:C43459"] = sample["SAMPLE_NAME"].get("SCIENTIFIC_NAME", "")
    annots["sep:00195"]["obi:0100026"]["data:1179"] = sample["SAMPLE_NAME"].get("TAXON_ID", "")

    # Iterate through the sample attributes
    if "SAMPLE_ATTRIBUTES" in sample and "SAMPLE_ATTRIBUTE" in sample["SAMPLE_ATTRIBUTES"]:
      attrs = sample["SAMPLE_ATTRIBUTES"]["SAMPLE_ATTRIBUTE"]
      if not isinstance(attrs, list):
        attrs = [attrs]

      for attr in attrs:
        # Skip tags with missing values.
        if attr["VALUE"].lower() == "missing":
            continue
        # Handle the cultivar
        elif attr["TAG"].lower() == "cultivar":
            annots["sep:00195"]["obi:0100026"]["GEMmaker:infraspecific_type"] = "cultivar"
            annots["sep:00195"]["obi:0100026"]["TAXRANK:0000045"] = attr["VALUE"]
        # Handle tags in the mapping file.
        elif attr["TAG"].lower() in tags:
            cvs = cv_map.loc[attr["TAG"].lower()]['CV_IDs'].split(',')
            code = "annots['%s'] = '%s'" % ("']['".join(cvs), attr["VALUE"])
            exec(code)
            print(code)
        else:
          sys.stderr.write("Unhandled sample attribute: \"%s\": \"%s\"\n" % (attr["TAG"], attr["VALUE"]))

    return annots







def get_accession_directory(meta_dir, accession):
    """
    Given an NCBI SRA accession number, generate a directory path for it.

    This function breaks each accession into 3-character string and
    creates a directory heirarchy based on the sequence of those sets of
    characters. If the directory does not exist it will be created.

    :param meta_dir:  The base directory where metadata files will be saved.
    :param accession: The NCBI SRA accession ID (experiment, run, sample, etc.)
    """
    dir = meta_dir + '/SRA_META/' + '/'.join([accession[i:i+3] for i in range(0, len(accession), 3)])
    pathlib.Path(dir).mkdir(parents=True, exist_ok=True)
    return dir






def create_amx(meta_dir, exp_ids, amx, cv_map):
    """
    Performs a looking for each run to get it's experiment and prints it.

    Prints the results to STDOUT

    :param mata_dir: The location where metafiles are stored.
    :param exp_ids: The list of SRA experiment IDs.
    """
    for exp_id in exp_ids:
        run = None
        experiment = None
        sample = None

        exp_dir = get_accession_directory(meta_dir, exp_id)
        exp_path = "%s/%s.ncbi.meta.json" % (exp_dir, exp_id)
        if os.path.exists(exp_path):
            with open(exp_path, "r") as exp_file:
                experiment = json.load(exp_file)

        if experiment is not None:
            sample_id = experiment['SAMPLE']['@accession']
            sample_dir = get_accession_directory(meta_dir, sample_id)
            sample_path = "%s/%s.ncbi.meta.json" % (sample_dir, sample_id)
            if os.path.exists(sample_path):
                with open(sample_path, "r") as sample_file:
                    sample = json.load(sample_file)

        # Convert the data to a JSON array of controlled vocabulary terms
        if experiment is not None and sample is not None:
            annots = map_annotations(experiment, sample, meta_dir, cv_map)
            print(annots)




if __name__ == "__main__":

    # Holds the mapping of SRA runs to experiments.
    run_to_exp = {}

    # Parse command-line arguments
    parser = argparse.ArgumentParser()
    parser.add_argument("--emx", help="The GEMmaker produced GEM file.", required=True)
    parser.add_argument("--amx", help="The file name for the output annotation matrix.", required=True)
    parser.add_argument("--meta_dir", help="The directory where the SRA_META folder is found", required=True)
    parser.add_argument("--cv_map", help="The path to the file that contains the mapping of NCBI SRA tags to controlled vocabulary terms. A default file is provided in the 'input/SRAannots2CV.txt' file of GEMmaker.")

    # Get the input arguments.
    args = parser.parse_args()
    meta_dir = args.meta_dir
    amx = args.amx

    # Read in the expression matrix and get the first row. These are the
    # experiment IDs.
    emx = pd.read_csv(args.emx, sep="\t", index_col=0, na_values = ['NA'])
    exp_ids = emx.columns

    # Read in the file that maps NCBI SRA tags to controlled vocabularies.
    cv_map = pd.read_csv(args.cv_map, sep="\t", comment='#')
    cv_map.index = cv_map['SRA_Tag']

    # Make the Metadata FAIR
    create_amx(meta_dir, exp_ids, amx, cv_map)
