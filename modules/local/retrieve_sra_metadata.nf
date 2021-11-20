/**
 * Retrieves metadata for all of the remote samples
 * and maps SRA runs to SRA experiments.
 */
process retrieve_sra_metadata {
  publishDir "${params.outdir}/reports", mode: params.publish_dir_mode
  container "systemsgenetics/gemmaker:2.0.0"

  input:
  file(sras)

  output:
  path("SRA_run2exp.tsv"), emit: SRR2SRX
  path("failed_runs.metadata.txt"), emit: FAILED_RUNS

  script:
  """
  >&2 echo "#TRACE n_remote_run_ids=`cat ${sras} | wc -l`"

  # Remove the 'path:' prefix. This was added to prevent
  # Nextflow from recoginzing the path and noticing the work
  # directory changed and trying to re-run this process even
  # if it succeeded.
  workdir=`echo ${params.workDirParent}/${params.workDirName}/GEMmaker | sed 's/path://'`

  retrieve_sra_metadata.py \
      --run_id_file ${sras} \
      --meta_dir "\$workdir" \
      --out_file SRA_run2exp.tsv \
       ${params.skip_samples ? "--skip_file ${file(params.skip_samples)}" : ""}
  """
}
