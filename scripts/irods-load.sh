#!/bin/bash
# Load the input files for an experiment from iRODS.

# parse command-line arguments
if [[ $# != 2 ]]; then
        echo "usage: $0 <remote-dir> <local-dir>"
        exit -1
fi

REMOTE_PATH="$1"
LOCAL_PATH="$2"

# load iRODS module
shopt -s expand_aliases

module load irods

# connect to iRODS server
iinit

# download input files
iget -Pr $REMOTE_PATH $LOCAL_PATH

# exit iRODS
iexit
