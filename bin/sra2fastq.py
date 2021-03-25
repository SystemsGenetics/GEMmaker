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
    failed_runs = {}

    # Now dump each SRA file
    for sra in sras:
        print("Dumping SRA {}".format(sra), file=sys.stdout)
        p = subprocess.Popen(["fastq-dump", "--split-files", sra],
                             stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        res = process_wait(p)

        # If the dump fails then lets mark this run as failed and move on.
        if (res['exit'] != 0):
            failed_runs[sra.replace('.sra', '')] = res['stderr']

    return(failed_runs)






if __name__ == "__main__":

    # Parse command-line arguments.
    parser = argparse.ArgumentParser()
    parser.add_argument("--sra_files", dest='sra_files', type=str, required=True,
                        help="List of input SRA files", nargs="+")
    parser.add_argument("--sample", dest='sample', type=str, required=True,
                        help="The sample name to which the SRA files belong")

    args = parser.parse_args()

    # Dump the samples. The names of those that failed get returned.
    failed_runs = dump_samples(args.sra_files)

    # Write any failed SRRs to a file
    f = open('{}.failed_runs.fastq-dump.txt'.format(args.sample), "w")
    f.write(json.dumps(failed_runs, indent=4))
    f.close()

    # Clean up any failed runs.
    if (len(failed_runs.keys()) > 0):
        print("Cleaning...", file=sys.stderr)
        fastq_files = glob.glob('*.fastq')
        for fastq_file in fastq_files:
            os.remove(fastq_file)

        # Create the file 'sample_failed' to trigger the workflow to skip
        # this sample.
        f = open('sample_failed', "w")
        f.close()

    # Return the exit code.
    sys.exit(0)
