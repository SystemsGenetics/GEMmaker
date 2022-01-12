/**
 * Performs KALLISTO alignemnt of fastq files
 */
process kallisto {
    tag { sample_id }
    publishDir "${params.outdir}/Samples/${sample_id}", mode: params.publish_dir_mode, pattern: params.publish_pattern_kallisto_ga

    conda (params.enable_conda ? "bioconda::kallisto=0.46.2" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/kallisto:0.46.2--h4f7b962_1"
    } else {
        container "quay.io/biocontainers/kallisto:0.46.2--h4f7b962_1"
    }

    input:
    tuple val(sample_id), path(fastq_files)
    path(kallisto_index)

    output:
    tuple val(sample_id), path("*.ga", type: "dir"), emit: GA_FILES
    tuple val(sample_id), path("*.kallisto.log"), emit: LOGS
    tuple val(sample_id), val(params.DONE_SENTINEL), emit: DONE_SIGNAL

    script:
    """
    echo "#TRACE sample_id=${sample_id}"
    echo "#TRACE fastq_lines=`cat *.fastq | wc -l`"

    # Convert the incoming FASTQ file list to an array
    fastq_files=(${fastq_files})

    # Kallisto will generate an exit code of 1 if no alignments are made. This
    # isn't an error and should be rewritten to 0 so that Nextflow doesn't end.
    trap 'if [[ \$? == 1 ]]; then echo OK; exit 0; fi' EXIT

    if [ \${#fastq_files[@]} == 2 ]; then
      echo "#TRACE Paired End Files Detected"
      kallisto quant \
        -i ${kallisto_index} \
        -b ${params.kallisto_bootstrap_samples} \
        -o ${sample_id}.Kallisto.ga \
        -t ${task.cpus} \
        \${fastq_files[0]} \
        \${fastq_files[1]} > ${sample_id}.kallisto.log 2>&1
    else
      echo "#TRACE Unpaired End Files Detected"
      kallisto quant \
        --single \
        -l 70 \
        -s .0000001 \
        -i ${kallisto_index} \
        -b ${params.kallisto_bootstrap_samples} \
        -o ${sample_id}.Kallisto.ga \
        -t ${task.cpus} \
        \${fastq_files[0]} > ${sample_id}.kallisto.log 2>&1
    fi
    """
}
