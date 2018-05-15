#!/usr/bin/env nextflow

/*
 * ------------------
 * SRA 2 GEV Pipeline
 * ------------------
 *
 * Authors:
 *  + John Hadish
 *  + Tyler Biggs
 *  + Stephen Ficklin
 *
 * Summary:
 *   A workflow for processing a large amount of SRA data...
 */


// Display the workflow title and show user parameters.
println """\
=================================
 S R A 2 G E V   P I P E L I N E
=================================

Parameters:
  + Remote SRA list path:        ${params.sra_list_path}
  + Local sample glob:           ${params.local_samples_path}
  + Genome reference path:       ${params.ref.path}
  + Reference genome prefix:     ${params.ref.prefix}
  + Trimmomatic clip path:       ${params.trimmomatic.clip_path}
  + Trimmomatic minimum length:  ${params.trimmomatic.MINLEN}
"""


/*
 * Local Sample Input.
 * This checks the folder that the user has given
 */

if (params.local_samples_path == 'none'){
  Channel
    .empty()
    .set { LOCAL_SAMPLES }
} else{
  Channel
    .fromFilePairs( params.local_samples_path, size: -1 )
    .set { LOCAL_SAMPLES }
}

// LOCAL_SAMPLES.subscribe{ println it }

/*
 * Remote SRA Input.
 */
Channel
  .from( file(params.sra_list_path).readLines() )
  // .filter { it }
  .set { REMOTE_SRAS }


/*
 * The fastq dump process downloads any needed remote fasta files to the
 * current working directory.
 */
process fastq_dump {

  module 'sratoolkit'
  publishDir "$sra", mode: 'link'
  time '24h'
  tag { sra }

  input:
    val sra from REMOTE_SRAS

  output:
    set val(sra), file("${sra}_?.fastq") into DOWNLOADED_SRAS

  """
    fastq-dump --split-files $sra
  """
}

// DOWNLOADED_SRAS.subscribe{ println it }

/*
 * Combine the remote and local samples into the same channel.
 */
COMBINED_SAMPLES = DOWNLOADED_SRAS.mix( LOCAL_SAMPLES )

// COMBINED_SAMPLES.subscribe{ println it }


/*
 * Performs a SRR to SRX converison:
 *
 * This first checks to see if the format is standard SRR,ERR,DRR
 * This takes the input SRR numbersd and converts them to SRX.
 * This is done by a python script that is stored in the "scripts" dir
 * The next step combines them
 */
process SRR_to_SRX {

  module 'python3'
  publishDir "$sra", mode: 'link'
  tag { sra }

  input:
    set val(sra), file(pass_files) from COMBINED_SAMPLES

  output:
    set stdout, file(pass_files) into SRX_GROUPS mode flatten

  """
  if [[ "$sra" == [SDE]RR* ]]; then
    python3 ${PWD}/scripts/retrieve_sample_metadata.py $sra

  else
    echo -n "SRX_$sra"
  fi
  """
}
// SRX_GROUPS.subscribe{ println it }

/*
 * This groups the channels based on srx numbers.
 */
SRX_GROUPS
  .groupTuple()
  .set { GROUPED_SRX }


// GROUPED_SRX.subscribe{ println it }

  /**
   *
   * This process merges the fastq files based on their SRX number.
   */
process SRR_combine{
  publishDir "$srx", mode: 'link'
  tag { srx }

  input:
    set val(srx), file(grouped) from GROUPED_SRX
  output:
    set val(srx), file("${srx}_?.fastq") into MERGED_SAMPLES

  // This command tests to see if ls produces a 0 or not by checking
  // its standard out. We do not use a "if [-e *foo]" becuase it gets
  // confused if there are more than one things returned by the wildcard
  """
    if ls *_1.fastq >/dev/null 2>&1; then
      cat *_1.fastq >> "${srx}_1.fastq"
    fi

    if ls *_2.fastq >/dev/null 2>&1; then
      cat *_2.fastq >> "${srx}_2.fastq"
    fi
  """

}



/*
 * Performs Trimmomatic on all fastq files.
 *
 * This process requires that the ILLUMINACLIP_PATH environment
 * variable be set in the trimmomatic module. This indicates
 * the path where the clipping files are stored.
 *
 */
 process trimmomatic {

   module "trimmomatic"
   publishDir "$srx", mode: 'link'
   tag { srx }

   input:
     set val(srx), file("${srx}_?.fastq") from MERGED_SAMPLES

   output:
     set val(srx), file("${srx}_?.trim.fastq"), file("${srx}_?s.trim.fastq") into TRIMMED_SAMPLES

   script:
       """
       if [ -e ${srx}_1.fastq ] && [ -e ${srx}_2.fastq ]; then
         java -Xmx512m org.usadellab.trimmomatic.Trimmomatic \
           PE \
           -threads 1 \
           -phred33 \
           ${srx}_1.fastq \
           ${srx}_2.fastq \
           ${srx}_1.trim.fastq \
           ${srx}_1s.trim.fastq \
           ${srx}_2.trim.fastq \
           ${srx}_2s.trim.fastq \
           ILLUMINACLIP:${params.trimmomatic.clip_path}/fasta_adapter.txt:2:40:15 \
           LEADING:3 \
           TRAILING:6 \
           SLIDINGWINDOW:4:15 \
           MINLEN:${params.trimmomatic.MINLEN}
       else
         # For ease of the next steps, rename the reverse file to the forward.
         # since these are non-paired it really shouldn't matter.
         if [ -e ${srx}_2.fastq]; then
           mv ${srx}_2.fastq ${srx}_1.fastq
         fi
         # Even though this is not paired-end, we need to create the 1s.trim.fastq
         # file as an empty file so that the rest of the workflow works
         touch ${srx}_1s.trim.fastq
         # Now run trimmomatic
         java -Xmx512m org.usadellab.trimmomatic.Trimmomatic \
           SE \
           -threads 1 \
           -phred33 \
           ${srx}_1.fastq \
           ${srx}_1.trim.fastq \
           ILLUMINACLIP:${params.trimmomatic.clip_path}/fasta_adapter.txt:2:40:15 \
           LEADING:3 \
           TRAILING:6 \
           SLIDINGWINDOW:4:15 \
           MINLEN:${params.trimmomatic.MINLEN}
       fi
       """
 }


/*
 * Performs hisat2 alignment of fastq files to a genome reference
 *
 * depends: trimmomatic
 */
process hisat2 {

  module 'hisat2'
  publishDir "$srx", mode: 'link'
  stageInMode "link"
  tag { srx }

  input:
   set val(srx), file("${srx}_?.trim.fastq"), file("${srx}_?s.trim.fastq") from TRIMMED_SAMPLES

  output:
   set val(srx), file("${srx}_vs_${params.ref.prefix}.sam") into INDEXED_SAMPLES

  script:
   """
     export HISAT2_INDEXES=${params.ref.path}
     if [ -e ${srx}_2.trim.fastq ]; then
       hisat2 \
         -x ${params.ref.prefix} \
         --no-spliced-alignment \
         -q \
         -1 ${srx}_1.trim.fastq \
         -2 ${srx}_2.trim.fastq \
         -U ${srx}_1s.trim.fastq,${srx}_2s.trim.fastq \
         -S ${srx}_vs_${params.ref.prefix}.sam \
         -t \
         -p 1 \
         --dta-cufflinks
     else
       hisat2 \
         -x ${params.ref.prefix} \
         --no-spliced-alignment \
         -q \
         -U ${srx}_1.trim.fastq \
         -S ${srx}_vs_${params.ref.prefix}.sam \
         -t \
         -p 1 \
         --dta-cufflinks
     fi
   """
}


/*
 * Sorts the SAM alignment file and coverts it to binary BAM
 *
 * depends: hisat2
 */
process samtools_sort {
  module 'samtools'
  publishDir "$srx", mode: 'link'
  stageInMode "link"
  tag { srx }

  input:
    set val(srx), file("${srx}_vs_${params.ref.prefix}.sam") from INDEXED_SAMPLES

  output:
    set val(srx), file("${srx}_vs_${params.ref.prefix}.bam") into SORTED_FOR_INDEX

  script:
    """
    samtools sort -o ${srx}_vs_${params.ref.prefix}.bam -O bam ${srx}_vs_${params.ref.prefix}.sam
    """
}


/*
 * Indexes the BAM alignment file
 *
 * depends: samtools_index
 */
process samtools_index {
  module 'samtools'
  publishDir "$srx", mode: 'link'
  stageInMode "link"
  tag { srx }

  input:
    set val(srx), file("${srx}_vs_${params.ref.prefix}.bam") from SORTED_FOR_INDEX

  output:
    set val(srx), file("${srx}_vs_${params.ref.prefix}.bam"), file("${srx}_vs_${params.ref.prefix}.bam.bai") into BAM_INDEXED_FOR_STRINGTIE

  script:
    """
    samtools index ${srx}_vs_${params.ref.prefix}.bam
    """
}



/**
 * Generates expression-level transcript abundance
 *
 * depends: samtools_index
 */
process stringtie {
  module 'stringtie'
  publishDir "$srx", mode: 'link'
  stageInMode "link"
  tag { srx }

  input:
    // We don't really need the .bai file, but we want to ensure
    // this process runs after the samtools_index step so we
    // require it as an input file.
    set val(srx), file("${srx}_vs_${params.ref.prefix}.bam"), file("${srx}_vs_${params.ref.prefix}.bam.bai") from BAM_INDEXED_FOR_STRINGTIE

  output:
    set val(srx), file("${srx}_vs_${params.ref.prefix}.ga") into STRINGTIE_GTF

  script:
    """
    stringtie \
    -v \
    -p 1 \
    -e \
    -o ${srx}_vs_${params.ref.prefix}.gtf \
    -G ${params.ref.path}/${params.ref.prefix}.gtf \
    -A ${srx}_vs_${params.ref.prefix}.ga \
    -l ${srx} ${srx}_vs_${params.ref.prefix}.bam
    """
}
/*
 * Generates the final FPKM file
 */
process fpkm {
  publishDir "$srx", mode: 'link'
  stageInMode "link"
  tag { srx }

  input:
    set val(srx), file("${srx}_vs_${params.ref.prefix}.ga") from STRINGTIE_GTF

  output:
    file "${srx}_vs_${params.ref.prefix}.fpkm" into FPKMS

  script:
    """
    ${PWD}/scripts/gtf2fpkm.sh ${srx} ${params.ref.prefix}
    """
}
