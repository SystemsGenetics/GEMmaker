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
for dir in ${directory}
do
  pwd
  echo Before Perl \$dir
  dir=`echo \$dir | perl -pi -e 's/[\\[,\\]]//g'`
  echo After Perl \$dir
  if [ -e \$dir ]; then
    echo \$dir
    pwd
    file_list=`ls -d \$dir/{*,}/{*,}`
    echo File List to delete: \$file_list

    for file in \$file_list
    do
      file=`echo \$file | perl -pi -e 's/[\\[,\\]]//g'`
      if [ -e \$file ]; then
        # Log some info about the file for debugging purposes
        echo "cleaning \$file"
        stat \$file
        # Get file info: size, access and modify times
        size=`stat --printf="%s" \$file`
        atime=`stat --printf="%X" \$file`
        mtime=`stat --printf="%Y" \$file`
        # Make the file size 0 and set as a sparse file
        > \$file
        truncate -s \$size \$file
        # Reset the timestamps on the file
        touch -a -d @\$atime \$file
        touch -m -d @\$mtime \$file
      fi
    done
    pwd
  fi
done
