/**
 * Downloads SRA files from NCBI using the SRA Toolkit.
 */
process download_runs {
    tag { sample_id }
    publishDir params.outdir, mode: params.publish_dir_mode, pattern: '*.failed_runs.download.txt', saveAs: { "Samples/${sample_id}/${it}" }
    container "systemsgenetics/gemmaker:2.0.0"

    input:
    tuple val(sample_id), val(run_ids), val(type)

    output:
    tuple val(sample_id), path("*.sra"), optional: true, emit: SRA_FILES
    tuple val(sample_id), path("sample_failed"), optional: true, emit: FAILED_SAMPLES
    tuple val(sample_id), path("*.failed_runs.download.txt"), emit: FAILED_RUNS

    script:
    """
    # Remove the 'path:' prefix. This was added to prevent
    # Nextflow from recoginzing the path and noticing the work
    # directory changed and trying to re-run this process even
    # if it succeeded.
    workdir=`echo ${params.workDirParent}/${params.workDirName}/GEMmaker | sed 's/path://'`

    echo "#TRACE sample_id=${sample_id}"
    echo "#TRACE n_remote_run_ids=${run_ids.tokenize(" ").size()}"
    echo "#TRACE n_spots=`retrieve_sra_spots.py \$workdir ${sample_id}`"

    retrieve_sra.py --sample ${sample_id} --run_ids ${run_ids} --akey \${ASPERA_KEY}
    """
}
