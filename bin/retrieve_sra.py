#!/usr/bin/env python3

"""
A Python script for retrieving a set of NCBI SRA run.

This script is meant as a replacement for the SRAToolkit
prefetch program, as it would sometimes fail without providing
proper errors. It requires Aspera and the SRAToolkit.

.. module:: GEMmaker
    :platform: UNIX, Linux
    :synopsis: This script receives a single input argument: a comma separated list of NCBI SRA run (SRR) IDs.
"""
import argparse
import sys
import subprocess
import re
import time
import os
import random
import requests
import shutil
import json

def process_wait(p):
    """
    A helper function for running a system command.

    :param p: A process that has already been opened.
    """

    print("{}".format(p.args), file=sys.stdout)

    # Wait for the process to finish to get the exit code,
    # STDOUT, and STDERR.  The STDOUT and STDERR need to be convereted
    # to a string.
    exit_code = p.wait()
    (stdout, stderr) = p.communicate()
    if (stdout):
        stdout = stdout.decode('utf-8')
    else:
        stdout = ''
    if (stderr):
        stderr = stderr.decode('utf-8')
    else:
        stderr = ''
    
    # Send stdout and stderr out to proper streams
    print(stdout, file=sys.stdout)
    print(stderr, file=sys.stderr)
   
    return { 'exit' : exit_code, 'stdout' : stdout, 'stderr' : stderr }

def get_sample_url(run_id):
    """
    Uses the srapath tool to get the URL of an SRA ID.

    :param run_id: the run ID.
    """

    print("Getting download paths for sample: {}".format(run_id))

    # First try if there is an aspera path.
    p = subprocess.Popen(["srapath", "--protocol", "fasp", "--json", run_id], stdout=subprocess.PIPE)
    res = process_wait(p)
    res = json.loads(res['stdout'])
    fasp_path = res['responses'][0]['remote'][0]['path']
    sra_size = res['responses'][0]['size']
 
    # If an https path was returned for Aspera then use that. 
    https_path = ""
    if (re.match('https://', fasp_path)): 
       https_path = fasp_path
       fasp_path = ""
    # If there is no aspera path then try HTTPs 
    elif (not fasp_path): 
       p = subprocess.Popen(["srapath", "--protocol", "https", "-P", run_id], stdout=subprocess.PIPE)
       res = process_wait(p)
       res = json.loads(res['stdout'])
       https_path = res['response'][0]['remote'][0]['path']
       sra_size = res['responses'][0]['size']

    urls = { 'https' : https_path, 'fasp' : fasp_path, 'size' : sra_size }
    return urls

def download_aspera(run_id, urls):
    """
    Downloads a SRA Run using Aspera.

    :param run_id:  the run ID.
    :param urls: a dictionary of urls for the run.
    """

    ec = 0
    print("Retrieving sample via aspera: {}".format(run_id))
    p = subprocess.Popen(["ascp", "-i", "$ASPERA_KEY","-k", "1", "-T", "-l", "1000m", urls['fasp'].replace('fasp://',''), "{}.sra".format(run_id)], stdout=subprocess.PIPE)
    res = process_wait(p)
    if (res['exit'] != 0):
        retry = False
    else:
        print("Aspera Failed. Exit code: {}. Trying https".format(exit_code), file=sys.stderr)
        ec = download_https(urls)
    return ec

def download_https(run_id, urls):
    """
    Downloads a SRA Run using HTTPS.

    :param run_id:  the run ID.
    :param urls: a dictionary of urls for the run.
    """

    ec = 0
    print("Retrieving sample via https: {}".format(run_id))
    try: 
        r = requests.get(urls['https'], verify=True, stream=True)
        r.raw.decode_content = True
        with open("{}.sra".format(run_id), 'wb') as sra_file:
            shutil.copyfileobj(r.raw, sra_file) 
    except requests.exceptions.RequestException as e:
        print(str(e), file=sys.stderr)
        ec = 1
    return ec 

def download_samples(run_ids):
    """
    Downloads a set of SRA Runs.

    :param run_ids: the list of run IDs.
    """

    # Set an exit code
    ec = 0

    # Now download each SRR.
    for run_id in run_ids:
        # Get the Aspera (fasp)and https URLs.
        # We will try Aspera first if we have a URL.
        urls = get_sample_url(run_id)
        if (urls['fasp']):
            ec = download_aspera(run_id, urls)
        else:
            ec = download_https(run_id, urls)
        if (ec != 0):
            print("Download failed.", file=sys.stderr)
            break
        
        if (sample_is_good(run_id, urls['size']) == False):
            print("Downloaded sample is missing or corrupted.", file=sys.stderr)
            ec = 1
            break

    return(ec)

def sample_is_good(run_id, size):
    """
    Checks if a sample is fully downloaded. 

    :param run_ids: the list of run IDs.
    """
    sra_file = "{}.sra".format(run_id)
    if (os.path.exists(sra_file)):
      statinfo = os.stat(sra_file)
      if (statinfo.st_size == size):
        print("Size {} == {}".format(statinfo.st_size, size), file=sys.stdout)
        return True
      else: 
        print("Size {} != {}".format(statinfo.st_size, size), file=sys.stdout)
    return False


if __name__ == "__main__":

    # Parse command-line arguments.
    parser = argparse.ArgumentParser()
    parser.add_argument("SRR_IDs", help="SRA run ID list provided as a single comma-separated string.")

    args = parser.parse_args()

    # Convert the SRR IDs into a list.
    run_ids = args.SRR_IDs.split(",")

    # Use this RE to make sure that each SRA run ID is correct.
    srr_re = re.compile("^[SED]RR\d+$")

    # Iterate through the run IDs and make sure they are all good.
    for run_id in run_ids:
        if not srr_re.match(run_id):
            raise ValueError("Improper SRA run ID: %s" % (run_id))

    # Download the samples:
    ec = download_samples(run_ids)

    # If the exit code is not zero then clean up so that
    # we don't waste space as Nextflow doesn't clean up
    # failed processes
    if (ec != 0):
        print("Cleaning after failed attempt.", file=sys.stderr)
        p = subprocess.Popen(["rm", "-rf", "./*"])

    # Return the exit code.
    sys.exit(ec)
