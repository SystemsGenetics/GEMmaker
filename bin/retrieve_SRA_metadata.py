#!/usr/bin/env python3

"""A Python script for retrieving metadata about NCBI SRA experiment runsself.

.. module:: GEMmaker
    :platform: UNIX, Linux
    :synopsis: This script recieves a single input argument: a file containing
       a list of run IDs (numbers with SRR, ERR or DRR prefixes) from NCBI's
       sequence read archive (SRA). It generates a variety of JSON meta files
       and outputs a tab-delimited file that maps run IDs to experiment IDs.
"""

import sys
import urllib
import xmltodict
import re
import pprint
import pandas as pd
import json

Script, File = sys.argv

def download_runs_meta(srr_file):
    """ Downloads the meta data for the runs contained the input file.

    :param srr_file:  the path of the file containing the run IDs. Each run
      ID must be on a new line in the file. Empty lines are ignored.
    :type srr_file: A string

    """
    # Holds the list of valid SRR Ids
    srr_ids = [];

    # Holds the list of experiment IDs and the their runs
    srx_ids = {};

    # The regular expression to ensure the SRA Run ID is correct
    sra_match = re.compile("^[SED]RR\d+$")

    # Open the SRR list file and make sure they all are good.
    fh = open(srr_file, "r")
    lines = fh.readlines();
    for line in lines:
      # Remove surrounding white space and skip empty lines.
      line = line.strip()
      if (not line):
          continue

      # If the SRR is valid then add it to our query list.
      if (sra_match.match(line)):
          srr_ids.append(line)
      else:
          raise ValueError('Improper SRR: "' + line + '"')

    # Open the file that we'll use to save the mapping of runs to experiments
    mapping_file = "SRR2SRX.txt"
    with open(mapping_file, 'w') as tabfile:

        # Now that we have our Ids we can make the bulk query. But there is a
        # limit of 10,000 at a time so we have to "page" our queries.
        page_size = 100;
        pages = int(len(srr_ids) / page_size + 1)
        for page in range(0, pages):

            url = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi'
            data = urllib.parse.urlencode({
              'retmod' : 'xml',
              'db' : 'sra',
              'retstart' : page * page_size,
              'retmax' : page_size,
              'id' : ','.join(srr_ids)
            }).encode()
            request = urllib.request.Request(url, data)
            response_obj = urllib.request.urlopen(request)
            response_xml = response_obj.read()
            response = xmltodict.parse(response_xml)

            # Get the list of experiments from the query.
            experiments = response["EXPERIMENT_PACKAGE_SET"]["EXPERIMENT_PACKAGE"]

            # If we only have one experiment then we have to convert it to an 
            # array for our for loop below.
            if (isinstance(experiments, list) == False):
              experiments = []
              experiments.append(response["EXPERIMENT_PACKAGE_SET"]["EXPERIMENT_PACKAGE"])

            # Now loop through the experiments to handle the metadata.
            for i in range(0, len(experiments)):
                experiment = experiments[i]
                run_id = srr_ids[page * page_size + i]
                exp_id = experiment["EXPERIMENT"]["IDENTIFIERS"]["PRIMARY_ID"];
                if exp_id in srx_ids:
                  srx_ids[exp_id].append(run_id)
                else:
                  srx_ids[exp_id] = [run_id];

                # Get the metadata for the sample.
                sample = experiment["SAMPLE"]

                # Get this run's metadata.
                # First, create the list of runs.
                run = [];
                runs = experiment["RUN_SET"]["RUN"]
                if isinstance(runs, list):
                  pass
                else:
                  runs = []
                  runs.append(experiment["RUN_SET"]["RUN"])
                # Next iterate through the runs until we find the right one.
                for item in runs:
                  if item["@accession"] == run_id:
                    run = item

                # Write out the experiment details in JSON metadata files.
                save_ncbi_meta(experiment, sample, run)

                # Now convert the data to a JSON array of controlled vocabulary terms
                save_gemmaker_meta(experiment, sample, run)

                # Now save the mapping information
                tabfile.write(run_id + "\t" + exp_id + "\n");

                # Write in CSV format to STDOUT for Nextflow to handle
                sys.stdout.write(run_id + "," + exp_id + "\n");


def save_ncbi_meta(experiment, sample, run):
    """ Creates .meta.json files containing run, experiment and sample info.

    :param experiment:  A dictionary containing the experiment metadata.
    :type experiment: dictionary

    :param sample: A dictionary containing the sample metadata.
    :type sample: dictionary

    :param run: A dictionary containing the run metadta.
    :type run: dictionaty.

    """
    # The query returns the experiment details. Write out
    # the metadata for the experiment.
    exp_id = experiment["EXPERIMENT"]["IDENTIFIERS"]["PRIMARY_ID"];
    with open(exp_id + '.ncbi.meta.json', 'w') as expfile:
      json.dump(experiment, expfile)

    # Next, write out the sample metadata file
    sample_id = sample['@accession']
    with open(sample_id + '.ncbi.meta.json', 'w') as samplefile:
      json.dump(sample, samplefile)

    # Finally, write out the run metadata file.
    run_id = run['@accession'];
    with open(run_id + '.ncbi.meta.json', 'w') as runfile:
      json.dump(run, runfile)


def save_gemmaker_meta(experiment, sample, run):
    """Writes a file using controlled vocabulary terms for meta data.

    This function creates both a JSON and tab delimited meta data file that
    maps meta data from NCBI to known controlled vocabulary terms. The purpose
    of this is to help ensure uniformity in meta data details between
    runs of GEMmaker and between different sample sets.

    :param experiment:  A dictionary containing the experiment metadata.
    :type experiment: dictionary

    :param sample:  A dictionary containing the sample metadata.
    :type sample: dictionary

    :param run:  A dictionary containing the run metatdata.
    :type: dictionary
    """
    # Save metadata using controlled vocabulary terms
    run_annots = pd.DataFrame()

    # Now that we have all the metadata loaded create the non-nested
    # annotation dictionary
    annots = {}

    # Run info
    annots['data:2091'] = run['@accession']

    # The Experiment
    annots['local:SRX_id'] = experiment['EXPERIMENT']['@accession']

    # Biological Sample
    annots['sep:00195'] = {}
    annots['sep:00195']['data:2091'] = sample['@accession']
    annots['sep:00195']['schema:title'] = sample['TITLE']
    annots['sep:00195']['schema:name'] = sample['@alias']
    annots['sep:00195']['obi:organism'] =  {}
    annots['sep:00195']['obi:organism']['rdfs:label'] =  sample['SAMPLE_NAME']['SCIENTIFIC_NAME']
    annots['sep:00195']['obi:organism']['NCIT:C43459'] =  sample['SAMPLE_NAME']['SCIENTIFIC_NAME']
    annots['sep:00195']['obi:organism']['data:1179'] =  sample['SAMPLE_NAME']['TAXON_ID']

    # Set defaults
    annots['sep:00195']['NCIT:C25150'] = ''
    annots['sep:00195']['NCIT:C16631'] = ''
    annots['sep:00195']['NCIT:C12801'] = ''
    annots['sep:00195']['NCIT:C43531'] = ''
    annots['NCIT:C25206'] = ''
    annots['EFO:0000721'] = ''
    annots['EFO:0000727'] = ''

    # Iterate through the sample attributes
    if 'SAMPLE_ATTRIBUTES' in sample and 'SAMPLE_ATTRIBUTE' in sample['SAMPLE_ATTRIBUTES']:
      attrs = sample['SAMPLE_ATTRIBUTES']['SAMPLE_ATTRIBUTE']
      for attr in attrs:

        # Add the cultivar
        if attr['TAG'] == 'cultivar':
          if attr['VALUE'] != 'missing':
            annots['sep:00195']['obi:organism']['local:infraspecific_type'] = 'cultivar'
            annots['sep:00195']['obi:organism']['TAXRANK:0000045'] = attr['VALUE']
          continue

        # Add the age
        if attr['TAG'] == 'age':
          if attr['VALUE'] != 'missing':
            annots['sep:00195']['NCIT:C25150'] = attr['VALUE']
          continue

        # Add the genotype
        if attr['TAG'] == 'Genotype' or attr['TAG'] == 'genotype':
          if attr['VALUE'] != 'missing':
            annots['sep:00195']['NCIT:C16631'] = attr['VALUE']
          continue

        # Add the tissue
        if attr['TAG'] == 'tissue':
          if attr['VALUE'] != 'missing':
            annots['sep:00195']['NCIT:C12801'] = attr['VALUE']
          continue

        # Add the developmental stage
        if attr['TAG'] == 'dev_stage':
          if attr['VALUE'] != 'missing':
            annots['sep:00195']['NCIT:C43531'] = attr['VALUE']
          continue

        # Add the temperature
        if attr['TAG'] == 'temp':
          if attr['VALUE'] != 'missing':
            annots['NCIT:C25206'] = attr['VALUE']
          continue

        # Add the time
        if attr['TAG'] == 'time':
          if attr['VALUE'] != 'missing':
            annots['EFO:0000721'] = attr['VALUE']
          continue

        # Add the treatment
        if attr['TAG'] == 'treatment':
          if attr['VALUE'] != 'missing':
            annots['EFO:0000727'] = attr['VALUE']
          continue

        sys.stderr.write("Unhandled sample attribute: '" + attr['TAG'] + "', value: " + attr['VALUE'] + "\n")

    # Save the heirarchical JSON metadata from GEMmaker
    jsonfilename = annots['local:SRX_id'] + '.GEMmaker.meta.json'
    with open(jsonfilename, 'w') as jsonfile:
      json.dump(annots, jsonfile)

    # Save a flattened CSV file
    tabfilename = annots['local:SRX_id'] + '.GEMmaker.meta.tab'
    with open(tabfilename, 'w') as tabfile:
      tabfile.write(
        annots["data:2091"]
        + "\t" + annots["local:SRX_id"]
        + "\t" + annots["sep:00195"]["data:2091"]
        + "\t" + annots["sep:00195"]["schema:title"]
        + "\t" + annots["sep:00195"]["schema:name"]
        + "\t" + annots["sep:00195"]["obi:organism"]["NCIT:C43459"]
        + "\t" + annots["sep:00195"]["NCIT:C25150"]
        + "\t" + annots["sep:00195"]["NCIT:C43531"]
        + "\t" + annots["sep:00195"]["NCIT:C16631"]
        + "\t" + annots["sep:00195"]["NCIT:C12801"]
        + "\t" + annots["NCIT:C25206"]
        + "\t" + annots["EFO:0000721"]
        + "\t" + annots['EFO:0000727']
        + "\n"
      )

if __name__ == "__main__":
    download_runs_meta(File)
