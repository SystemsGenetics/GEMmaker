#!/usr/bin/env python3

"""
A Python script for retrieving a set of NCBI SRA run.

This script is meant to gracefully handle transfer errors that can occur
when downloading SRR files (e.g. timeouts, incomplete downloads, etc.).  If an
error does occur then depending on the prefetch exit code, the fetch is retried.
The script will always return 0 if the SRRs are all downloaded successfully or
1 if not. If any of the SRRs fail, the directory is emptied so as not to overrun
storage on successive runs.

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

def download_samples(run_ids):
    """
    Downloads a set of SRA Runs.

    :param run_ids: the list of run IDs.
    """

    # Set an exit code
    ec = 0

    # Now download each SRR.
    for run_id in run_ids:
        num_retries = 0;
        retry = True
        max_retries = 5;
        while retry == True:
            print("Retrieving sample: {}".format(run_id))
            
            # Run Prefetch with support for Aspera. This expects that the
            # ascp program is in the PATH and that there is an Environment
            # variable naemd $ASPERA_KEY that has the path to the SSH key.
            p = subprocess.Popen(["prefetch", "-v", "--max-size", "50G", "--output-directory", ".", "--ascp-path", "`which ascp`\"|$ASPERA_KEY\"", "--ascp-options", "-k 1 -T -l 1000m", run_id], stdout=subprocess.PIPE)

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

            # Send stdout and stderr out to calling program
            print(stdout, file=sys.stdout)
            print(stderr, file=sys.stderr)

            # Check to see if there were any common problems.
            if (exit_code == 0):
                retry = False
            else:
                print("Failed. Exit code: {}".format(exit_code), file=sys.stderr)

                # Exit code #3: transfer incomplete while reading file within
                # network system module
                if (exit_code == 3):
                    print("Transfer incomplete.  sleeping for a bit and then trying again...", file=sys.stderr)
                    p = subprocess.Popen(["rm", "-rf", run_id])
                    # Sleep for 10 minutes to give things time to "cool off"
                    time.sleep(600)
                # If we've encountered an exit code that we're not familiar
                # with then exit. We can then add new code here to address
                # other codes.
                else:
                    retry = False;
                    ec = 1

                # If retry is set to true then only allow it for up to
                # max_retries
                if (retry == True):
                    num_retries = num_retries + 1
                    if num_retries >= max_retries:
                        print("Exceeded max retries.", file=sys.stderr)
                        retry = False
                        ec = 1
                    else:
                        print("Retry number {}".format(num_retries), file=sys.stderr)

        # Sometimes prefetch will download the SRA into a set of files and
        # sometimes into a directory. It's probably a setting somewhere but we
        # need consistency. So move the .sra files into the working directory.
        if (ec == 0 and os.path.exists("{}".format(run_id))):
            if (os.path.exists("./{}/{}.sra".format(run_id, run_id))):
                subprocess.Popen(["mv", "{}/{}.sra".format(run_id, run_id), "."])
            if (os.path.exists("./{}/{}_1.sra".format(run_id, run_id))):
                subprocess.Popen(["mv", "{}/{}_1.sra".format(run_id, run_id), "."])
                subprocess.Popen(["mv", "{}/{}_2.sra".format(run_id, run_id), "."])
    return(ec)

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
