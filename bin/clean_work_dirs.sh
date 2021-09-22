#!/bin/bash
# This script is meant for cleaning all files out of a directory in a Nextflow
# Work directory. This will investigate the directory up to 2 deep, and replace
# all files with a sparse file so that they will not take up space.
#
# The file cleaning process empties the file, converts it to a sparse file so it
# has an acutal size of zero but appears as the original size, the access
# and modify times are kept the same.
#
# This is based off of the clean_work_files.sh template
# Creator: John Hadish
directory="$1"

for dir in ${directory}; do
  if [ -e $dir ]; then
    echo "Cleaning: $dir"
    files=`find $dir -type  f `

    echo "Files to delete: $files"
    clean_work_files.sh "$files" "null"
  fi
done
