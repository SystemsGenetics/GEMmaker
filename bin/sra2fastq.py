#!/usr/bin/env python3

"""
A Python script for dumping FASTQ files from an SRA

.. module:: GEMmaker
    :platform: UNIX, Linux
    :synopsis: This script receives a single input argument: a comma separated list of NCBI SRA run (SRR) IDs.
"""
import argparse
import sys
import subprocess
import os
import glob
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

def dump_samples(sras):
    """
    dumps a SRA files to FASTQ files

    :param sras: the list of SRA file
    """

    # Store the list of failed SRA files.
    failed = {}

    # Now dump each SRA file
    for sra in sras:
        print("Dumping SRA {}".format(sra), file=sys.stdout)
        p = subprocess.Popen(["fastq-dump", "--split-files", sra],
                             stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        res = process_wait(p)

        # If the dump fails then lets mark this run as failed and move on.
        if (res['exit'] != 0):
            failed[sra] = res['stderr']
            print("Cleaning up after failed attempt.", file=sys.stderr)

            # Remove the FASTQ files
            fastq_list = glob.glob('{}_*.fastq'.format(sra.replace('sra','')))
            for fastq in fastq_list:
                os.remove(fastq)

    return(failed)

def redownload_srr(srr):
  """Re-downloads an SRA file.

  :param srr: The name of the SRR file to download
  :type srr: string

  This is needed in the event that an SRR file is corrupted and we want to try
  again.  This will force the file to be downloaded into the current directory.
  """
  p = subprocess.Popen(["../../../bin/retrieve_sra.py", srr],
                       stdout=subprocess.PIPE, stderr=subprocess.PIPE)
  return process_wait(p)


if __name__ == "__main__":

    # Parse command-line arguments.
    parser = argparse.ArgumentParser()
    parser.add_argument("sra_files", help="list of input SRA files", nargs="+")

    args = parser.parse_args()

    # Dump the samples. The names of those that failed get returned.
    failed = dump_samples(args.sra_files)

    # Write any failed SRRs to a file
    f = open('failed_runs.txt', "w")
    f.write(json.dumps(failed, indent=4))
    f.close()

    # Return the exit code.
    sys.exit(0)
