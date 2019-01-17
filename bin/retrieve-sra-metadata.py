#!/usr/bin/env python3

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



def download_runs_meta(run_ids, page_size=100):
    """
    Downloads the metadata for the runs contained the input file.

    :param run_ids: the list of run IDs
    """
    experiments = []

    # query the run IDs in "pages" because the NCBI endpoint is fragile
    for idx in range(0, len(run_ids), page_size):
        sys.stderr.write("Downloading run IDs %6d - %6d of %6d...\n" % (idx, idx + page_size - 1, len(run_ids)))

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

        # parse the XML response
        response_obj = urllib.request.urlopen(request)
        response_xml = response_obj.read()
        response = xmltodict.parse(response_xml)

        # get the list of experiments from the query
        page_experiments = response["EXPERIMENT_PACKAGE_SET"]["EXPERIMENT_PACKAGE"]

        # append page experiments to experiments list
        if isinstance(page_experiments, list):
            experiments += page_experiments
        else:
            experiments.append(page_experiments)

    # process the metadata from each experiment
    for run_id, experiment in zip(run_ids, experiments):
        # Get the experiment ID
        exp_id = experiment["EXPERIMENT"]["IDENTIFIERS"]["PRIMARY_ID"]

        # Get the metadata for the sample
        sample = experiment["SAMPLE"]

        # Get the list of runs
        runs = experiment["RUN_SET"]["RUN"]
        if not isinstance(runs, list):
            runs = [runs]

        # Get the metadata for this run
        run = []
        for item in runs:
            if item["@accession"] == run_id:
                run = item

        # Write out the experiment details in JSON metadata files
        save_ncbi_meta(experiment, sample, run)

        # Convert the data to a JSON array of controlled vocabulary terms
        save_gemmaker_meta(experiment, sample, run)

        # Write in CSV format to stdout
        sys.stdout.write("%s,%s\n" % (run_id, exp_id))



def save_ncbi_meta(experiment, sample, run):
    """
    Creates .meta.json files containing run, experiment and sample info.

    :param experiment: A dictionary containing the experiment metadata.
    :param sample: A dictionary containing the sample metadata.
    :param run: A dictionary containing the run metadta.
    """
    # Write out the metadata for the experiment
    exp_id = experiment["EXPERIMENT"]["IDENTIFIERS"]["PRIMARY_ID"]
    expfile = open("%s.ncbi.meta.json" % (exp_id), "w")
    json.dump(experiment, expfile)

    # Write out the sample metadata file
    sample_id = sample["@accession"]
    samplefile = open("%s.ncbi.meta.json" % (sample_id), "w")
    json.dump(sample, samplefile)

    # Write out the run metadata file
    run_id = run["@accession"]
    runfile = open("%s.ncbi.meta.json" % (run_id), "w")
    json.dump(run, runfile)



def save_gemmaker_meta(experiment, sample, run):
    """
    Writes a file using controlled vocabulary terms for metadata.

    This function creates both a JSON and tab delimited metadata file that
    maps metadata from NCBI to known controlled vocabulary terms. The purpose
    of this is to help ensure uniformity in metadata details between
    runs of GEMmaker and between different sample sets.

    :param experiment: A dictionary containing the experiment metadata.
    :param sample: A dictionary containing the sample metadata.
    :param run: A dictionary containing the run metatdata.
    """
    # Now that we have all the metadata loaded create the non-nested
    # annotation dictionary
    annots = {}

    # Run info
    annots["data:2091"] = run["@accession"]

    # The Experiment
    annots["local:SRX_id"] = experiment["EXPERIMENT"]["@accession"]

    # Biological Sample
    annots["sep:00195"] = {}
    annots["sep:00195"]["data:2091"] = sample.get("@accession", "")
    annots["sep:00195"]["schema:title"] = sample.get("TITLE", "")
    annots["sep:00195"]["schema:name"] = sample.get("@alias", "")
    annots["sep:00195"]["obi:organism"] = {}
    annots["sep:00195"]["obi:organism"]["rdfs:label"] = sample["SAMPLE_NAME"].get("SCIENTIFIC_NAME", "")
    annots["sep:00195"]["obi:organism"]["NCIT:C43459"] = sample["SAMPLE_NAME"].get("SCIENTIFIC_NAME", "")
    annots["sep:00195"]["obi:organism"]["data:1179"] = sample["SAMPLE_NAME"].get("TAXON_ID", "")

    # Set defaults
    annots["sep:00195"]["NCIT:C25150"] = ""
    annots["sep:00195"]["NCIT:C16631"] = ""
    annots["sep:00195"]["NCIT:C12801"] = ""
    annots["sep:00195"]["NCIT:C43531"] = ""
    annots["NCIT:C25206"] = ""
    annots["EFO:0000721"] = ""
    annots["EFO:0000727"] = ""

    # Iterate through the sample attributes
    if "SAMPLE_ATTRIBUTES" in sample and "SAMPLE_ATTRIBUTE" in sample["SAMPLE_ATTRIBUTES"]:
      attrs = sample["SAMPLE_ATTRIBUTES"]["SAMPLE_ATTRIBUTE"]
      if not isinstance(attrs, list):
        attrs = [attrs]

      for attr in attrs:
        # Add the cultivar
        if attr["TAG"] == "cultivar":
          if attr["VALUE"] != "missing":
            annots["sep:00195"]["obi:organism"]["local:infraspecific_type"] = "cultivar"
            annots["sep:00195"]["obi:organism"]["TAXRANK:0000045"] = attr["VALUE"]

        # Add the age
        elif attr["TAG"] == "age":
          if attr["VALUE"] != "missing":
            annots["sep:00195"]["NCIT:C25150"] = attr["VALUE"]

        # Add the genotype
        elif attr["TAG"] == "Genotype" or attr["TAG"] == "genotype":
          if attr["VALUE"] != "missing":
            annots["sep:00195"]["NCIT:C16631"] = attr["VALUE"]

        # Add the tissue
        elif attr["TAG"] == "tissue":
          if attr["VALUE"] != "missing":
            annots["sep:00195"]["NCIT:C12801"] = attr["VALUE"]

        # Add the developmental stage
        elif attr["TAG"] == "dev_stage":
          if attr["VALUE"] != "missing":
            annots["sep:00195"]["NCIT:C43531"] = attr["VALUE"]

        # Add the temperature
        elif attr["TAG"] == "temp":
          if attr["VALUE"] != "missing":
            annots["NCIT:C25206"] = attr["VALUE"]

        # Add the time
        elif attr["TAG"] == "time":
          if attr["VALUE"] != "missing":
            annots["EFO:0000721"] = attr["VALUE"]

        # Add the treatment
        elif attr["TAG"] == "treatment":
          if attr["VALUE"] != "missing":
            annots["EFO:0000727"] = attr["VALUE"]

        else:
          sys.stderr.write("Unhandled sample attribute: \"%s\": \"%s\"\n" % (attr["TAG"], attr["VALUE"]))

    # Save the heirarchical JSON metadata from GEMmaker
    jsonfilename = "%s.GEMmaker.meta.json" % (annots["local:SRX_id"])
    jsonfile = open(jsonfilename, "w")
    json.dump(annots, jsonfile)

    # Save a flattened CSV file
    tabfilename = "%s.GEMmaker.meta.tab" % (annots["local:SRX_id"])
    tabfile = codecs.open(tabfilename, "w", "utf-8")

    fields = [
        annots["data:2091"],
        annots["local:SRX_id"],
        annots["sep:00195"]["data:2091"],
        annots["sep:00195"]["schema:title"],
        annots["sep:00195"]["schema:name"],
        annots["sep:00195"]["obi:organism"]["NCIT:C43459"],
        annots["sep:00195"]["NCIT:C25150"],
        annots["sep:00195"]["NCIT:C43531"],
        annots["sep:00195"]["NCIT:C16631"],
        annots["sep:00195"]["NCIT:C12801"],
        annots["NCIT:C25206"],
        annots["EFO:0000721"],
        annots["EFO:0000727"]
    ]

    tabfile.write("\t".join(fields) + "\n")



if __name__ == "__main__":
    # parse command-line arguments
    parser = argparse.ArgumentParser()
    parser.add_argument("SRR_ID_FILE", help="SRA run ID file")
    parser.add_argument("--page-size", type=int, default=100, help="number of SRA run IDs to query at a time", dest="PAGE_SIZE")

    args = parser.parse_args()

    # load SRA run IDs
    srr_file = open(args.SRR_ID_FILE, "r")
    run_ids = [line.strip() for line in srr_file]
    run_ids = [run_id for run_id in run_ids if run_id]

    # make sure that each SRA run ID is correct
    sra_pattern = re.compile("^[SED]RR\d+$")

    for run_id in run_ids:
      if not sra_pattern.match(run_id):
          raise ValueError("Improper SRA run ID: %s" % (run_id))

    # download metadata for each SRA run ID
    download_runs_meta(run_ids, page_size=args.PAGE_SIZE)
