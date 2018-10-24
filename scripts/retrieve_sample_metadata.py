import sys
import urllib
import xmltodict

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
    sys.stdout.write(
        response
        ["EXPERIMENT_PACKAGE_SET"]
        ["EXPERIMENT_PACKAGE"]
        ["EXPERIMENT"]
        ["IDENTIFIERS"]
        ["PRIMARY_ID"])


# -----------------------------------------------------------------------------
# Code below is run if this file is called as a script.
# -----------------------------------------------------------------------------
if __name__ == "__main__":
    download_json(SRR)
