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
    clean_work_files.sh "${files.join(" ")}"
    """
}
