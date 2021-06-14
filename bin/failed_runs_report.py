#!/usr/bin/env python3

"""
A Python script for generated an HTML file

.. module:: GEMmaker
    :platform: UNIX, Linux
    :synopsis: This script receives a single input argument: a comma separated list of NCBI SRA run (SRR) IDs.
"""
import argparse
import glob
import json
import sys
import jinja2
import os





def get_failed_metadata():
    items = []
    if os.path.exists('failed_runs.metadata.txt'):
        with open('failed_runs.metadata.txt') as f:
            data = json.load(f)
            for item in data.keys():
                items.append({'run_id': item, 'reason': data[item]})
    return(items)






def get_failed_downloads():
    items = []
    file_list = glob.glob('*.failed_runs.download.txt')
    for file in file_list:
        sample_id = file.replace('.failed_runs.download.txt', '')
        with open(file) as f:
            data = json.load(f)
            for item in data.keys():
                items.append({'run_id': item,
                              'sample_id': sample_id,
                              'reason': data[item]})
    return(items)





def get_failed_dumps():
    items = []
    file_list = glob.glob('*.failed_runs.fastq-dump.txt')
    for file in file_list:
        sample_id = file.replace('.failed_runs.fastq-dump.txt', '')
        with open(file) as f:
            data = json.load(f)
            for item in data.keys():
                items.append({'run_id': item,
                              'sample_id': sample_id,
                              'reason': data[item]})
    return(items)


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("--template", dest='template', type=str, required=True,
                        help="Path to the report tempalte file")
    args = parser.parse_args()

    failed_metadata = get_failed_metadata()
    failed_downloads = get_failed_downloads()
    failed_dumps = get_failed_dumps()


    subs = jinja2.Environment(loader=jinja2.FileSystemLoader('./')).\
      get_template(args.template).\
      render(title="GEMmaker Failed Runs Report",
             failed_metadata = failed_metadata,
             failed_downloads = failed_downloads,
             failed_dumps = failed_dumps
             )

    # lets write the substitution to a file
    with open('failed_SRA_run_report.html','w') as f: f.write(subs)

    # Return the exit code.
    sys.exit(0)
