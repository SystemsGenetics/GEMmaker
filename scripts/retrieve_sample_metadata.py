import sys
import urllib
import xmltodict
import json

query = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=sra&term='
fetch = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=sra&id='

script, SRR = sys.argv

def download_json (line):
    '''
    This function downlads the data for a SRR and returns its SRX.
    '''
    # Get the current SRA ID.
    sra_id = line.strip()

    # First perform the entrez query with the SRA database.
    response_obj = urllib.request.urlopen(query + sra_id + '[Accession]')
    response_xml = response_obj.read()
    response = xmltodict.parse(response_xml)
    query_id = response['eSearchResult']['IdList']['Id']

    # Next get the query results. We are only querying a single SRA
    # record at a time.
    response_obj = urllib.request.urlopen(fetch + query_id)
    response_xml = response_obj.read()
    response = xmltodict.parse(response_xml)

    # This returns all the SRR numbers associated with this SRX
    SRX = response["EXPERIMENT_PACKAGE_SET"]["EXPERIMENT_PACKAGE"]["EXPERIMENT"]["IDENTIFIERS"]["PRIMARY_ID"]
    sys.stdout.write(SRX)

    # The query returns the experiment details. Write out
    # the metadata for the experiment.
    exp = response["EXPERIMENT_PACKAGE_SET"]["EXPERIMENT_PACKAGE"]
    with open(SRX + '.ncbi.meta.json', 'w') as expfile:
      json.dump(exp, expfile)

    # Next, write out the sample meta data file
    sample = response["EXPERIMENT_PACKAGE_SET"]["EXPERIMENT_PACKAGE"]["SAMPLE"]
    sample_id = sample['@accession']
    with open(sample_id + '.ncbi.meta.json', 'w') as samplefile:
      json.dump(sample, samplefile)

    # Next, find the run details and write that meta data file
    runs = response["EXPERIMENT_PACKAGE_SET"]["EXPERIMENT_PACKAGE"]["RUN_SET"]["RUN"]
    if isinstance(runs, list):
      pass
    else:
      runs = []
      runs.append(response["EXPERIMENT_PACKAGE_SET"]["EXPERIMENT_PACKAGE"]["RUN_SET"]["RUN"])

    for run in runs:
      if run["@accession"] == SRR:
        with open(SRR + '.ncbi.meta.json', 'w') as runfile:
          json.dump(run, runfile)


if __name__ == "__main__":
    download_json(SRR)
