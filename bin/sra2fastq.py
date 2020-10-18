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
import re
import time
import os
import random
import requests
import shutil

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

    # Set an exit code
    ec = 0

    # Now dump each SRA
    for sra in sras:
        print("Dumping SRA {}".format(sra), file=sys.stdout)
        p = subprocess.Popen(["fastq-dump", "--split-files", sra], stdout=subprocess.PIPE)
        res = process_wait(p)
        ec = res['exit']
        if (ec != 0):
            ec = 1
            print("Dump to FASTQ failed.", file=sys.stderr)
            break

    return(ec)


if __name__ == "__main__":

    # Parse command-line arguments.
    parser = argparse.ArgumentParser()
    parser.add_argument("sra_files", help="list of input SRA files", nargs="+")

    args = parser.parse_args()

    # Download the samples:
    ec = dump_samples(args.sra_files)

    # If the exit code is not zero then clean up so that
    # we don't waste space as Nextflow doesn't clean up
    # failed processes
    if (ec != 0):
        print("Cleaning after failed attempt.", file=sys.stderr)
        p = subprocess.Popen(["rm", "-rf", "./*"])

    # Return the exit code.
    sys.exit(ec)
