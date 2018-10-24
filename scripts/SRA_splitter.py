#!/usr/bin/env python

"""

This is a helper script for splitting a very large file
of NCBI SRAs into smaller chunk files for use with GEMmaker.

Usage:
.. code-block:: sh
    python SRA_splitter.py [file] [chunk-size]

The [file] argument should be the name of the file with the
list of SRA IDs. It should be a two-column tab-delimited file
where the fist column is the SRA experiment ID and the second
is the SRA run ID.

This script will then split the file into files with at most
[chunk-size] lines per file. Because GEMmaker merges all
runs for an experiment into a single sample, this script ensures
that the runs of an experiment are not split between two files.

"""

import argparse
import glob
import os
import pandas as pd
import sys
from textwrap import dedent
import logging

# -----------------------------------------------------------------------------
# Log configuration.
# -----------------------------------------------------------------------------
logging.basicConfig(level=logging.INFO, stream=sys.stdout)

# -----------------------------------------------------------------------------
# Set up command line argument parser.
# -----------------------------------------------------------------------------
parser = argparse.ArgumentParser(description=__doc__)

parser.add_argument('--sra-file', dest='SRA_file', action='store',
                    required=True, help=dedent("""\
    Name of the file with the list of SRA IDs. It should be a 
    two-column tab-delimited file where the fist column is the 
    SRA experiment ID and the second is the SRA run ID."""))

parser.add_argument('--chunk-size', dest='chunk_size', action='store',
                    type=int, required=True, help=dedent("""\
    Split the file into files with at most [chunk-size] lines 
    per file.."""))


# -----------------------------------------------------------------------------
# Separate main() function so that the code can be imported and tested.
# -----------------------------------------------------------------------------
def main(sra_file: str, chunk_size: int):
    """Split the SRA IDs found within the given SRA file into files of
    the given chunk size.

    :param sra_file: File which contains all the SRA IDs of interest.
    :param chunk_size: The size of SRA ID files to be produced.
    :returns: Writes new files of chunk size to disk.

    """

    def chunk_filename(chunk):
        """Helper function to build a file name by chunk value."""
        return "SRA_IDs" + f"{chunk:2.0}" + ".txt"

    # Remove existing SRA_IDs.*.txt files.
    for f in glob.glob("SRA_IDs.*.txt"):
        os.remove(f)

    # Load the SRA file as a pandas DataFrame.
    logging.info(f"Reading file {sra_file}.")

    samples = pd.read_table(sra_file, names=["exp_id", "run_id"])

    # Set the initial chunk number and current size.
    chunk_num = 1
    curr_size = 0

    # Create a file handle -- file handles are created and closed
    # separate from the logic of the loops below.
    file_handle = open(chunk_filename(chunk_num), "r")

    # Iterate through unique experiment IDs.
    for index, exp in sorted(set(samples["exp_id"])):

        # Fetch the run IDs for the given experiment ID.
        runs = samples.loc[samples["exp_id"] == exp, "run_id"]

        # Throw error if there are too many run IDs for chunk size.
        if len(runs) > chunk_size:
            logging.error(dedent(f"""Experiment {exp} has {len(runs)} \
                runs, which cannot be satisfied by chunk size of \
                {chunk_size}."""))
            file_handle.close()
            sys.exit(-1)

        while curr_size > chunk_size:

            # Determine whether the current chunk file has enough
            # lines remaining for the given run IDs.
            if curr_size + len(runs) > chunk_size:
                # Close the current file handle and log the number of runs.
                file_handle.close()

                logging.info(dedent(f"""Wrote {curr_size} runs \
                    to {chunk_filename(chunk_num)}."""))

                # Move on to the next chunk file.
                chunk_num += 1
                curr_size = 0
                file_handle = open(chunk_filename(chunk_num), "r")

            # Write a line per run to the current file handle, then increment
            # the count by that amount.
            file_handle.writelines(runs)
            curr_size += len(runs)

            logging.debug(f"Wrote {len(runs)} lines to {chunk_filename}.")

    # Close the last chunk file.
    logging.debug(f"Wrote {len(runs)} lines to {chunk_filename}.")
    file_handle.close()


# -----------------------------------------------------------------------------
# Code below is run if this file is called as a script.
# -----------------------------------------------------------------------------
if __name__ == "__main__":
    # Read in the input arguments and call the main() function.
    args = parser.parse_args()

    logging.info(dedent(f"""\
        {'=' * 80}
        {__name__} called from the command line with the following arguments:
        {vars(args)}"""))

    main(sra_file=args.SRA_file, chunk_size=args.chunk_size)
