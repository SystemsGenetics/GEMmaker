#!/usr/bin/env python

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

    stdout = stdout.decode('utf-8') if stdout else ''
    stderr = stderr.decode('utf-8') if stderr else ''

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

    # Initialize paths.
    fasp_path = ""
    https_path = ""

    # First try if there is an aspera path.
    p = subprocess.Popen(
        ["srapath", "--protocol", "fasp", "--json", run_id],
        stderr=subprocess.PIPE,
        stdout=subprocess.PIPE)

    res = process_wait(p)

    if res['exit'] == 0:
        res = json.loads(res['stdout'])
        fasp_path = res['responses'][0]['remote'][0]['path']
        sra_size = res['responses'][0]['size']
    else:
        print('Failed to fetch sample url via fasp, trying with http...', file=sys.stderr)

    # If an https path was returned for Aspera then use that.
    if (re.match('https://', fasp_path)):
        https_path = fasp_path
        fasp_path = ""

    # If there is no aspera path then try HTTPs
    elif (not fasp_path):
        p = subprocess.Popen(
            ["srapath", "--protocol", "https", "--json", run_id],
            stderr=subprocess.PIPE,
            stdout=subprocess.PIPE)

        res = process_wait(p)

        if res['exit'] == 0:
            res = json.loads(res['stdout'])
            https_path = res['responses'][0]['remote'][0]['path']
            sra_size = res['responses'][0]['size']
        else:
            print('Failed to fetch sample url via https.', file=sys.stderr)

    return { 'fasp' : fasp_path, 'https' : https_path, 'size' : sra_size }



def download_aspera(run_id, urls):
    """
    Downloads a SRA Run using Aspera.

    :param run_id:  the run ID.
    :param urls: a dictionary of urls for the run.
    """
    print("Retrieving sample via aspera: {}".format(run_id))

    p = subprocess.Popen(
        ["ascp", "-i", "$ASPERA_KEY","-k", "1", "-T", "-l", "1000m", urls['fasp'].replace('fasp://',''), "{}.sra".format(run_id)],
        stderr=subprocess.PIPE,
        stdout=subprocess.PIPE)

    res = process_wait(p)

    if (res['exit'] != 0):
        print("Aspera Failed: ".format(res['stderr']), file=sys.stderr)
        if 'https' in urls.keys():
            print("Trying https", file=sys.stderr)
            res = download_https(run_id, urls)

    return res



def download_https(run_id, urls):
    """
    Downloads a SRA Run using HTTPS.

    :param run_id:  the run ID.
    :param urls: a dictionary of urls for the run.
    """
    print("Retrieving sample via https: {}".format(run_id))

    res = {'exit': 0}

    try:
        r = requests.get(urls['https'], verify=True, stream=True)
        r.raw.decode_content = True

        with open("{}.sra".format(run_id), 'wb') as sra_file:
            shutil.copyfileobj(r.raw, sra_file)

    except requests.exceptions.RequestException as e:
        print(str(e), file=sys.stderr)
        res['exit'] = 1
        res['stderr'] = str(e)

    return res



def download_samples(run_ids):
    """
    Downloads a set of SRA Runs.

    :param run_ids: the list of run IDs.
    """
    failed_runs = {}

    # Now download each SRR.
    for run_id in run_ids:
        # Get the Aspera (fasp)and https URLs.
        # We will try Aspera first if we have a URL.
        urls = get_sample_url(run_id)

        if (urls['fasp']):
            res = download_aspera(run_id, urls)
        elif (urls['https']):
            res = download_https(run_id, urls)
        else:
            message = "Failed to fetch sample url."
            failed_run[run_id] = message
            print(message, file=sys.stderr)
            break

        if (res['exit'] != 0):
            message = "Download failed: {}".format(res['stderr'])
            failed_run[run_id] = message
            print(message, file=sys.stderr)
            break

        if (sample_is_good(run_id, urls['size']) == False):
            message = "Downloaded sample is missing or corrupted."
            failed_run[run_id] = message
            print(message, file=sys.stderr)
            break

    return failed_runs



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
    parser.add_argument("--run_ids", dest='run_ids', type=str, required=True,
                        help="List of input SRA files", nargs="+")
    parser.add_argument("--sample", dest='sample', type=str, required=True,
                        help="The sample name to which the SRA files belong")
    args = parser.parse_args()

    # Use this RE to make sure that each SRA run ID is correct.
    srr_re = re.compile("^[SED]RR\d+$")

    # Iterate through the run IDs and make sure they are all good.
    for run_id in args.run_ids:
        if not srr_re.match(run_id):
            raise ValueError("Improper SRA run ID: %s" % (run_id))

    # Download the samples:
    failed_runs = download_samples(args.run_ids)

    # Write any failed SRRs to a file
    f = open('{}.failed_runs.download.txt'.format(args.sample), "w")
    f.write(json.dumps(failed_runs, indent=4))
    f.close()

    # Clean up any failed runs.
    if (len(failed_runs.keys()) > 0):
        print("Cleaning....", file=sys.stderr)
        sra_files = glob.glob('*.sra')
        for sra_file in sra_files:
            os.remove(sra_file)

        # Create the file 'sample_failed' to trigger the workflow to skip
        # this sample.
        f = open('sample_failed', "w")
        f.close()

    # Return the exit code.
    sys.exit(0)
