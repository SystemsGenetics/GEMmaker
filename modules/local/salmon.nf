/**
 * Performs SALMON alignemnt of fastq files
 */
process salmon {
    tag { sample_id }
    publishDir "${params.outdir}/Samples/${sample_id}", mode: params.publish_dir_mode, pattern: params.publish_pattern_salmon_ga

    conda (params.enable_conda ? 'bioconda::salmon=1.5.2' : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/salmon:1.5.2--h84f40af_0"
    } else {
        container "quay.io/biocontainers/salmon:1.5.2--h84f40af_0"
    }

    input:
    tuple val(sample_id), path(fastq_files)
    path(salmon_index)

    output:
    tuple val(sample_id), path("*.ga", type: "dir"), emit: GA_FILES
    tuple val(sample_id), path("*.salmon.log"), emit: LOGS
    tuple val(sample_id), val(params.DONE_SENTINEL), emit: DONE_SIGNAL

    script:
    """
    echo "#TRACE sample_id=${sample_id}"
    echo "#TRACE fastq_lines=`cat *.fastq | wc -l`"

    # convert the incoming FASTQ file list to an array
    fastq_files=(${fastq_files})

    if [ \${#fastq_files[@]} == 2 ]; then
      salmon quant \
        -i ${salmon_index} \
        -l A \
        -1 \${fastq_files[0]} \
        -2 \${fastq_files[1]} \
        -p ${task.cpus} \
        -o ${sample_id}.Salmon.ga \
        --minAssignedFrags 1 > ${sample_id}.salmon.log 2>&1
    else
      salmon quant \
        -i ${salmon_index} \
        -l A \
        -r \${fastq_files[0]} \
        -p ${task.cpus} \
        -o ${sample_id}.Salmon.ga \
        --minAssignedFrags 1 > ${sample_id}.salmon.log 2>&1
    fi

    # Copy these files for MultiQC reporting
    cp ${sample_id}.Salmon.ga/aux_info/meta_info.json ${sample_id}-meta_info.json
    cp ${sample_id}.Salmon.ga/libParams/flenDist.txt ${sample_id}-flenDist.txt
    """
}
