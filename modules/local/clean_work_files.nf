// Import generic module functions
include { saveFiles } from './functions'

params.options = [:]

/**
 * Clean intermediate files.
 */
process clean_work_files {
    tag { sample_id }
    label "local"

    input:
    tuple val(sample_id), val(files)

    output:
    val(1), emit: IS_CLEAN

    script:
    """
    files="${files.join(" ")}"

    for file in \${files}; do
      if [ -e \$file ]; then
        # Log some info about the file for debugging purposes
        echo "Cleaning: \$file"
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
    """
}
