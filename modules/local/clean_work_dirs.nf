// Import generic module functions
include { saveFiles } from './functions'

params.options = [:]

/**
 * Clean intermediate directories.
 */
process clean_work_dirs {
    tag { sample_id }
    label "local"

    input:
    tuple val(sample_id), val(directory)

    script:
    """
    clean_work_dirs.sh "${directory}"
    """
}
