#!/usr/bin/env python

"""A Python script for creating gene expression matrix (GEM) files


.. module:: GEMmaker
  :platform: Unix, Windows
  :synopsis: The Gene expression matrix (GEM) file can be created after a
    successful run of GEMmaker.  It can be used to create GEM files containing
    either FPKM or TPM values.  It can also be used to combine the results
    from multiple GEMmaker directories into a single GEM file.

.. moduleauthor:: John Hadish & Stephen Ficklin

"""

import argparse
import logging
import pandas as pd
import glob
import os
import re
from textwrap import dedent

# -----------------------------------------------------------------------------
# Log configuration.
# -----------------------------------------------------------------------------
logging.basicConfig(level=logging.INFO, filemode="a",
                    filename="create_GEM.log")

# -----------------------------------------------------------------------------
# Specify the arguments that are allowed by this script
# -----------------------------------------------------------------------------
parser = argparse.ArgumentParser(description=dedent("""\
    The Gene expression matrix (GEM) file can be created after a
    successful run of GEMmaker.  It can be used to create GEM
    files containing either FPKM or TPM values.  It can also be 
    used to combine the results from multiple GEMmaker directories
    into a single GEM file."""))

parser.add_argument('--sources', dest='path', action='store', required=True,
                    nargs='*', help=dedent(""" One or more GEMmaker directory
    paths where results are stored. FPKM or TPM files are found in Sample_* or 
    [SDR]RX directories within the GEMmaker directory."""))

parser.add_argument('--prefix', dest='prefix', action='store', required=True,
                    help=dedent("""A prefix to be added to the GEM.txt file 
                    when writing the final matrix."""))

parser.add_argument('--type', dest='data_type', action='store', required=True,
                    choices=['TPM', 'FPKM'], help=dedent("""The type of count
                    values to include in the  GEM (either "TPM" or "FPKM"."""))


# -----------------------------------------------------------------------------
# Separate main() function so that the code can be imported and tested.
# -----------------------------------------------------------------------------
def main(prefix, data_type, path_list):
    """Combine FPKM or TPM values from one or more GEMmaker runs into
    a single GEM file.

    :param prefix: (str) The prefix of the .GEM file to be output.
    :param data_type: (str) Either "TPM" or "FPKM".
    :param path_list: (list(str)) A list of paths to examine for data_type.

    :returns: Constructs and saves to disk a single GEM file. Returns the
        path to which the file was saved.

    """

    # Set the name of the output GEM file.
    ematrix_name = prefix + ".GEM." + data_type + '.txt'

    # Iterate through the GEMmaker directories to find the FPKM and TPM files.
    result_files = []

    for source_dir in path_list:

        # Build the source directory glob.
        sra_glob = os.path.join(source_dir, "[SED]RX*", "*",
                                data_type.lower())

        logging.info(dedent(f"""Finding {data_type} files in {source_dir}.
            with {sra_glob}"""))

        # Check for NCBI SRA files.
        for filename in glob.iglob(sra_glob, recursive=True):
            result_files.append(filename)

        # Build a glob for non-SRA files.
        non_sra_glob = os.path.join(source_dir, "Sample_*",
                                    "*", data_type.lower())

        logging.info(dedent(f"""Finding {data_type} files in {source_dir}.
            with {non_sra_glob}"""))

        # Check for local non SRA files
        for filename in glob.iglob(non_sra_glob, recursive=True):
            result_files.append(filename)

    logging.info(dedent(f"""Found {len(result_files)} sample files."""))

    # Initialize the expression matrix by reading in the
    ematrix = pd.DataFrame({'gene': []})

    for result in result_files:
        # Get the sample name from the file name.
        file_basename = os.path.basename(result)
        sample_name = file_basename.split('_vs_')[0]

        # Remove the Sample_ from local file names.
        sample_name = re.sub(r'^Sample_', '', str(sample_name))

        logging.info(dedent(f"""Adding results for sample: {sample_name}"""))

        df = pd.read_csv(result, header=None, sep='\t',
                         names=["gene", sample_name])

        # Stringtie was creating duplicate entries for some genes because of
        # how it handles genes that have non-overlapping splice variants. So,
        # we need to remove the duplicates but keep the first occurrence.
        df.drop_duplicates(['gene'], keep='first', inplace=True)

        # Now merge in this sample to the expression matrix.
        ematrix = ematrix.merge(df.iloc[:, [0, 1]], on='gene', how='outer')

    # Set the gene names as the data frame indexes.
    ematrix = ematrix.set_index('gene', drop=True)

    # Write out our expression matrix.
    logging.info(dedent(f"""Writing {ematrix_name}."""))

    ematrix.to_csv(ematrix_name, sep='\t', na_rep="NA", index_label=False)


if __name__ == "__main__":
    # Read in the input arguments and call the main() function.
    args = parser.parse_args()
    main(path_list=args.path, prefix=args.prefix, data_type=args.data_type)
