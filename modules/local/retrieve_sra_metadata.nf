// Import generic module functions
include { saveFiles } from './functions'

params.options = [:]

/**
 * Retrieves metadata for all of the remote samples
 * and maps SRA runs to SRA experiments.
 */
process retrieve_sra_metadata {
  publishDir params.outdir, mode: params.publish_dir_mode, pattern: "failed_runs.metadata.txt"
  container "systemsgenetics/gemmaker:2.0.0"

  output:
  stdout emit: SRR2SRX
  path("failed_runs.metadata.txt"), emit: FAILED_RUNS

  script:
  """
  >&2 echo "#TRACE n_remote_run_ids=`cat ${file(params.sras)} | wc -l`"

  retrieve_sra_metadata.py \
      --run_id_file ${file(params.sras)} \
      --meta_dir ${workflow.workDir}/GEMmaker \
       ${params.skip_samples ? "--skip_file ${file(params.skip_samples)}" : ""}
  """
}
