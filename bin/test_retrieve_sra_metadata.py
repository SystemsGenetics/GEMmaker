import sys
from retrieve_sra_metadata import download_runs_meta

def test_download_runs_meta():

    # TEST #1: missing SRR.
    # The first item in this list does not exist in the NCBI SRA. To prevent
    # The request from failing we add in a second SRR that ensures we get
    # a response but we only have one run returned.
    test_set1 = ['SRR2927685', 'SRR4042625']
    download_runs_meta(test_set1)
