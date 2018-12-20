#!/bin/bash
# Save the output files from a experiment to iRODS.

# parse command-line arguments
if [[ $# != 2 ]]; then
        echo "usage: $0 <local-dir> <remote-dir>"
        exit -1
fi

LOCAL_PATH="$(realpath $1)"
REMOTE_PATH="$2"

# load iRODS module
shopt -s expand_aliases

module load irods

# connect to iRODS server
iinit

# upload output files
iput -bPr $LOCAL_PATH $REMOTE_PATH

# exit iRODS
iexit
