//
// This file holds several functions specific to the workflow/GEMmaker.nf in the systgemsgenetics/gemmaker pipeline
//
import java.io.File

class WorkflowGEMmaker {

    //
    // Check and validate parameters
    //
    public static void initialise(workflow, params, log) {

        // Create the directories used for running batches
        File gemmaker_dir = new File("${workflow.workDir}/GEMmaker")
        if (gemmaker_dir.isEmpty()) {
            gemmaker_dir.mkdir()
        }
        File stage_dir = new File("${workflow.workDir}/GEMmaker/stage")
        if (stage_dir.isEmpty()) {
            stage_dir.mkdir()
        }
        File process_dir = new File("${workflow.workDir}/GEMmaker/process")
        if (process_dir.isEmpty()) {
            process_dir.mkdir()
        }
        File done_dir = new File("${workflow.workDir}/GEMmaker/done")
        if (done_dir.isEmpty()) {
            done_dir.mkdir()
        }

        // Make sure that the user hasn't changed the quantification
        // tool on a resumed run.
        File method_lock_file = new File("${workflow.workDir}/GEMmaker/method")
        if (method_lock_file.isEmpty()) {
            method_lock_file << params.pipeline
        }
        else {
            def reader = method_lock_file.newReader()
            def active_method = reader.readLine()
            reader.close()
            if (!active_method.equals(params.pipeline)) {
                error "Error: previously, GEMmaker was set to run using the '${active_method}' tool, but it looks as though the configuration has changed to use the '${params.pipeline}' tool. GEMmaker only supports use of one tool at a time. If you would like to change the quantification tool please re-run GEMmaker in a new directory or remove the `work` and `results` directories prior to restarting GEMmaker to clear out unwanted results."
            }
        }
    }
}
