/**
 * This process merges the fastq files based on their sample_id number.
 */
process fastq_merge {
    tag { sample_id }
    container "systemsgenetics/gemmaker:2.0.0"

    input:
    tuple val(sample_id), path(fastq_files)

    output:
    tuple val(sample_id), path("${sample_id}_?.fastq"), emit: FASTQ_FILES
    tuple val(sample_id), val(params.DONE_SENTINEL), emit: DONE_SIGNAL

    script:
    """
    echo "#TRACE sample_id=${sample_id}"
    echo "#TRACE fastq_lines=`cat *.fastq | wc -l`"

    # First, concatenate all of the set 1 files
    files1=`ls *_1.fastq | grep -v ${sample_id} | sort`
    for file in $files1; do
       echo "Concatenate file: ${file} to ${sample_id}_1.fastq"
       cat $file >> "${sample_id}_1.fastq"
    done
    echo "Done with ${sample_id}_1.fastq"

    # Next, concatenate all of the set 2 files
    files2=`ls *_2.fastq | grep -v ${sample_id} | sort`
    for file in $files2; do
      echo "Concatenate file: ${file} to ${sample_id}_2.fastq"
      cat $file >> "${sample_id}_2.fastq"
    done
    echo "Done with ${sample_id}_2.fastq"

    # If there is a FASTQ sample with _2 suffix but no _1  then rename
    # to _1 so that count software will work
    if [ -e ${sample_id}_2.fastq ] && [ ! -e ${sample_id}_1.fastq ]; then
      mv ${sample_id}_2.fastq ${sample_id}_1.fastq
    fi
    """
}
