#!/usr/bin/env nextflow
/*
========================================================================================
    systemsgenetics/gemmaker
========================================================================================
    Github : https://github.com/systemsgenetics/gemmaker
    Website: https://gemmaker.readthedocs.io/en/latest/
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl = 2

/*
========================================================================================
    VALIDATE & PRINT PARAMETER SUMMARY
========================================================================================
*/

WorkflowMain.initialise(workflow, params, log)

/*
========================================================================================
    NAMED WORKFLOW FOR PIPELINE
========================================================================================
*/

include { GEMmaker } from './workflows/GEMmaker'

//
// WORKFLOW: Run main systemsgenetics/gemmaker analysis pipeline
//
workflow SYSTEMSGENETICS_GEMMAKER {
    GEMmaker()
}

/*
========================================================================================
    RUN ALL WORKFLOWS
========================================================================================
*/

//
// WORKFLOW: Execute a single named workflow for the pipeline
// See: https://github.com/nf-core/rnaseq/issues/619
//
workflow {
    GEMmaker()
}

/*
========================================================================================
    THE END
========================================================================================
*/
