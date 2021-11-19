// Import generic module functions
include { saveFiles } from './functions'

params.options = [:]

/**
 * Extracts FASTQ files from downloaded SRA files.
 */
process fastq_dump {
    tag { sample_id }
    publishDir params.outdir, mode: params.publish_dir_mode, pattern: params.publish_pattern_fastq_dump, saveAs: { "Samples/${sample_id}/${it}" }
    publishDir params.outdir, mode: params.publish_dir_mode, pattern: '*.failed_runs.fastq-dump.txt', saveAs: { "Samples/${sample_id}/${it}" }
    container "systemsgenetics/gemmaker:2.0.0"

    input:
    tuple val(sample_id), path(sra_files)

    output:
    tuple val(sample_id), path("*.fastq"), optional: true, emit: FASTQ_FILES
    tuple val(sample_id), path("sample_failed"), optional: true, emit: FAILED_SAMPLES
    tuple val(sample_id), path("*.failed_runs.fastq-dump.txt"), emit: FAILED_RUNS
    tuple val(sample_id), val(params.DONE_SENTINEL), emit: DONE_SIGNAL

    script:
    """
    echo "#TRACE sample_id=${sample_id}"
    echo "#TRACE sra_bytes=`stat -Lc '%s' *.sra | awk '{sum += \$1} END {print sum}'`"

    sra2fastq.py --sample ${sample_id} --sra_files ${sra_files}
    """
}
