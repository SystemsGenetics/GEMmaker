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

    script:
    """
    clean_work_files.sh "${files.join(" ")}"
    """
}
