#!/usr/bin/env python

"""
A Python script for retrieving metadata about NCBI SRA experiment runs.

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





def download_runs_meta(run_ids, meta_dir, page_size=100):
    """
    Downloads the metadata for the runs contained the input file.

    :param run_ids: the list of run IDs
    """
    experiments = []
    found_runs = []
    failed_runs = {}

    # query the run IDs in "pages" because the NCBI endpoint is fragile
    page = 0
    for idx in range(0, len(run_ids), page_size):
        sys.stderr.write("Retrieving metadata for run IDs %6d - %6d of %6d...\n" % (idx, idx + page_size - 1, len(run_ids)))

        # get the next page of IDs
        ids = run_ids[idx : (idx + page_size)]

        # query the next page of metadata
        url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
        data = urllib.parse.urlencode({
            "retmod": "xml",
            "db": "sra",
            "id": ",".join(ids)
        }).encode()
        request = urllib.request.Request(url, data)

        sys.stderr.write("Fetching IDs: %s.  " % (",".join(ids)))

        # parse the XML response
        try:
            response_obj = urllib.request.urlopen(request)
        except urllib.error.URLError as e:
            for id in ids:
                failed_runs[id] = e.reason
            sys.stderr.write("ERROR Retrieving SRA Metadata: %s\n" % (e.reason))
            continue

        response_xml = response_obj.read().decode(response_obj.headers.get_content_charset())
        response = xmltodict.parse(response_xml)

        # Write out the XML for debugging purposes should something fail
        outxml = "ncbi-sra-meta-%d.xml" % (page)
        xmlfile = open(outxml, "w", encoding='utf-8')
        xmlfile.write(response_xml)
        xmlfile.close()

        # make sure the XML is formatted correctly.
        if (not response["EXPERIMENT_PACKAGE_SET"]):
            raise Exception('XML returned from NCBI is missing EXPERIMENT_PACKAGE_SET.')
        if (not response["EXPERIMENT_PACKAGE_SET"]["EXPERIMENT_PACKAGE"]):
            raise Exception('XML returned from NCBI is missing EXPERIMENT_PACKAGE.')

        # get the list of experiments from the query
        page_experiments = response["EXPERIMENT_PACKAGE_SET"]["EXPERIMENT_PACKAGE"]

        # Make sure we have an array of experiments for looping over,
        # even if we only have one returned.
        if not isinstance(page_experiments, list):
            page_experiments = [page_experiments]

        experiments += page_experiments

        page = page + 1

    # remove duplicate experiments
    experiments = {exp["EXPERIMENT"]["IDENTIFIERS"]["PRIMARY_ID"]: exp for exp in experiments}.values()

    # process the metadata from each experiment
    n_runs_found = 0
    for experiment in experiments:

        # Get the experiment ID
        exp_id = experiment["EXPERIMENT"]["IDENTIFIERS"]["PRIMARY_ID"]

        # Get the metadata for the sample
        sample = experiment["SAMPLE"]

        # Get the list of runs
        runs = experiment["RUN_SET"]["RUN"]
        if not isinstance(runs, list):
            runs = [runs]

        # If the run belongs to an ID we passed into the SRA lookup
        # then we want to act on it.
        for run in runs:
            run_id = run["@accession"]
            if run_id in run_ids:
                n_runs_found = n_runs_found + 1
                found_runs.append(run_id)

                # Write out the experiment details in JSON metadata files
                save_ncbi_meta(experiment, sample, run, meta_dir)

            else:
                message = 'Notice: the run, %s, is part of experiment %s but was not included in the input SRA_IDS.txt file. Do you want to include this run?' % (run_id, exp_id)
                failed_runs[run_id] = message
                sys.stderr.write(message)

    if n_runs_found != len(run_ids):
        sys.stderr.write('Notice: could not retrieve metadata for the all runs: %d != %d' % (n_runs_found, len(run_ids)))
        bad_ids = set(run_ids).difference(set(found_runs))
        for bad_id in bad_ids:
            failed_runs[bad_id] = 'Metadata was not returned by NCBI for this run.'

    sys.stderr.write("Metadata for %d runs retrieved\n" % (n_runs_found))
    return(failed_runs)





def save_ncbi_meta(experiment, sample, run, meta_dir):
    """
    Creates .meta.json files containing run, experiment and sample info.

    :param experiment: A dictionary containing the experiment metadata.
    :param sample: A dictionary containing the sample metadata.
    :param run: A dictionary containing the run metadta.
    """
    # Write out the metadata for the experiment if it doesn't exist.
    exp_id = experiment["EXPERIMENT"]["IDENTIFIERS"]["PRIMARY_ID"]
    exp_dir = get_accession_directory(meta_dir, exp_id)
    exp_path = "%s/%s.ncbi.meta.json" % (exp_dir, exp_id)
    if not os.path.exists(exp_path):
        expfile = open(exp_path, "w")
        json.dump(experiment, expfile)
        expfile.close()

    # Write out the sample metadata file if it doesn't exist
    sample_id = sample["@accession"]
    sample_dir = get_accession_directory(meta_dir, sample_id)
    sample_path = "%s/%s.ncbi.meta.json" % (sample_dir, sample_id)
    if not os.path.exists(sample_path):
        samplefile = open(sample_path, "w")
        json.dump(sample, samplefile)
        samplefile.close()

    # Write out the run metadata file if it doesn't exist.
    run_id = run["@accession"]
    run_dir = get_accession_directory(meta_dir, run_id)
    run_path = "%s/%s.ncbi.meta.json" % (run_dir, run_id)
    if not os.path.exists(run_path):
        runfile = open(run_path, "w")
        json.dump(run, runfile)
        runfile.close()





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





def find_downloaded_runs(meta_dir, run_ids):
    """
    Checks if runs have already been downloaded.

    :param mata_dir: The location where metafiles are stored.
    :param run_ids: The list of SRA run IDs.

    :return: A list of run IDs that do not have metadata.
    """
    missing_runs = []

    for run_id in run_ids:
        run_dir = get_accession_directory(meta_dir, run_id)
        run_path = "%s/%s.ncbi.meta.json" % (run_dir, run_id)
        if not os.path.exists(run_path):
            missing_runs.append(run_id)

    return missing_runs





def map_run_to_exp(meta_dir, run_ids):
    """
    Performs a looking for each run to get it's experiment and prints it.

    Prints the results to STDOUT

    :param mata_dir: The location where metafiles are stored.
    :param run_ids: The list of SRA run IDs.
    """

    for run_id in run_ids:
        run_dir = get_accession_directory(meta_dir, run_id)
        run_path = "%s/%s.ncbi.meta.json" % (run_dir, run_id)
        if os.path.exists(run_path):
            with open(run_path, "r") as run_file:
                run = json.load(run_file)
                exp_id = run['EXPERIMENT_REF']['@accession']
                # Write in CSV format to stdout
                sys.stdout.write("%s,%s\n" % (run_id, exp_id))





if __name__ == "__main__":

    # Holds the mapping of SRA runs to experiments.
    run_to_exp = {}

    # Parse command-line arguments
    parser = argparse.ArgumentParser()
    parser.add_argument("--run_id_file", help="A file containing the SRA run IDs.", required=True)
    parser.add_argument("--skip_file", help="A file containing the SRA run IDs to skip", required=False)
    parser.add_argument("--meta_dir", help="The directory were meta data should be stored", required=True)
    parser.add_argument("--page-size", help="number of SRA run IDs to query at a time", type=int, default=100, dest="PAGE_SIZE")

    args = parser.parse_args()

    # Get the meta output dir
    meta_dir = args.meta_dir

    # Load SRA run IDs
    srr_file = open(args.run_id_file, "r")
    run_ids = [line.strip() for line in srr_file]
    run_ids = [run_id for run_id in run_ids if run_id]
    srr_file.close()

    # Load the file for skipping and remove any IDs to skip.
    if args.skip_file and os.path.exists(args.skip_file):
        skip_file = open(args.skip_file, "r")
        skip_ids = [line.strip() for line in skip_file]
        skip_ids = [skip_id for skip_id in skip_ids if skip_id]
        skip_file.close()
        run_ids = list(set(run_ids).difference(set(skip_ids)))

    # Make sure that each SRA run ID is correct
    sra_pattern = re.compile("^[SED]RR\d+$")

    for run_id in run_ids:
        if not sra_pattern.match(run_id):
            raise ValueError("Improper SRA run ID: %s" % (run_id))

    # Find runs whose metadata has already been retrieved.
    missing_runs = find_downloaded_runs(meta_dir, run_ids)

    if len(run_ids) - len(missing_runs) > 0:
        sys.stderr.write("Found %d SRA run file(s) already retreived. Skipping retreival of these.\n" % (len(run_ids) - len(missing_runs)))

    # Download metadata for each SRA run ID that is not already present
    failed_runs = download_runs_meta(missing_runs, meta_dir, args.PAGE_SIZE)

    # Saved failed runs into a file.
    f = open('failed_runs.metadata.txt', 'w')
    f.write(json.dumps(failed_runs, indent=4))
    f.close()

    # Now iterate through the metadata and and map runs to experiments.
    map_run_to_exp(meta_dir, run_ids)

    # Return the exit code.
    sys.exit(0)
