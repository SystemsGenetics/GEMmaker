
// Import generic module functions
include { saveFiles } from './functions'

params.options = [:]


/**
 * Move a new sample into the process directory when
 * a previous sample is completed.
 */
process next_sample {
    tag { sample_id }
    label "local"
    cache false
    maxForks 1

    input:
    val(sample_id)

    exec:
    // Move the completed file into the done folder.
    sample_file = file("${workflow.workDir}/GEMmaker/process/${sample_id}.sample.csv")
    sample_file.moveTo("${workflow.workDir}/GEMmaker/done")

    // Move the next sample file into the processing directory.
    staged_files = file("${workflow.workDir}/GEMmaker/stage/*")
    if (staged_files.size() > 0) {
        staged_files.first().moveTo("${workflow.workDir}/GEMmaker/process")
    }

    // Write the "done" file if there are no more samples to process.
    else {
        done_file = file("${workflow.workDir}/GEMmaker/process/DONE")
        done_file << ""
    }
}
