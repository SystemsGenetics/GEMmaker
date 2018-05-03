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
 */

if (params.local_samples_path == 'none'){
  Channel
    .empty()
    .set { LOCAL_SAMPLES}
} else{
  Channel
    .fromFilePairs( params.local_samples_path, size: -1 )
  // .i fEmpty()
    .set { LOCAL_SAMPLES }
}

/*
 * Remote SRA Input.
 */
Channel
  .from( file(params.sra_list_path).readLines() )
  .filter { it }
  .set { REMOTE_SRAS }


/*
 * The fastq dump process downloads any needed remote fasta files to the
 * current working directory.
 */
process fastq_dump {

  module 'sratoolkit'
  publishDir "$sra", mode: 'link'
  time '24h'

  input:
    val sra from REMOTE_SRAS

  output:
    set val(sra), file("${sra}_?.fastq") into DOWNLOADED_SRAS mode flatten

  """
    fastq-dump --split-files $sra
  """
}


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
    set val(sra), file(pass_files) from DOWNLOADED_SRAS

  output:
    set stdout, file(pass_files) into SRX_GROUPS

  """
  if [[ "$sra" == [SDE]RR* ]]; then
    python3 ${PWD}/scripts/retrieve_sample_metadata.py $sra
  else
    echo -n "SRX_$sra"
  fi
  """
}


/*
 * This groups the channels based on srx numbers.
 */
SRX_GROUPS
  .groupTuple()
  .set { GROUPED_SRX }


/*
 * Combine the remote and local samples into the same channel.
 */
COMBINED_SAMPLES = GROUPED_SRX.mix( LOCAL_SAMPLES )

// COMBINED_SAMPLES.subscribe{ println it }


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
   publishDir "$sra", mode: 'link'
   tag { sra }

   input:
     set val(sra), file("${sra}_?.fastq") from COMBINED_SAMPLES

   output:
     set val(sra), file("${sra}_?.trim.fastq"), file("${sra}_?s.trim.fastq") into TRIMMED_SAMPLES

   script:
       """
       if [ -e ${sra}_1.fastq ] && [ -e ${sra}_2.fastq ]; then
         java -Xmx512m org.usadellab.trimmomatic.Trimmomatic \
           PE \
           -threads 1 \
           -phred33 \
           ${sra}_1.fastq \
           ${sra}_2.fastq \
           ${sra}_1.trim.fastq \
           ${sra}_1s.trim.fastq \
           ${sra}_2.trim.fastq \
           ${sra}_2s.trim.fastq \
           ILLUMINACLIP:${params.trimmomatic.clip_path}/fasta_adapter.txt:2:40:15 \
           LEADING:3 \
           TRAILING:6 \
           SLIDINGWINDOW:4:15 \
           MINLEN:${params.trimmomatic.MINLEN}
       else
         # For ease of the next steps, rename the reverse file to the forward.
         # since these are non-paired it really shouldn't matter.
         if [ -e ${sra}_2.fastq]; then
           mv ${sra}_2.fastq ${sra}_1.fastq
         fi
         # Even though this is not paired-end, we need to create the 1s.trim.fastq
         # file as an empty file so that the rest of the workflow works
         touch ${sra}_1s.trim.fastq
         # Now run trimmomatic
         java -Xmx512m org.usadellab.trimmomatic.Trimmomatic \
           SE \
           -threads 1 \
           -phred33 \
           ${sra}_1.fastq \
           ${sra}_1.trim.fastq \
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
  publishDir "$sra", mode: 'link'
  stageInMode "link"
  tag { sra }

  input:
   set val(sra), file("${sra}_?.trim.fastq"), file("${sra}_?s.trim.fastq") from TRIMMED_SAMPLES

  output:
   set val(sra), file("${sra}_vs_${params.ref.prefix}.sam") into INDEXED_SAMPLES

  script:
   """
     export HISAT2_INDEXES=${params.ref.path}
     if [ -e ${sra}_2.trim.fastq ]; then
       hisat2 \
         -x ${params.ref.prefix} \
         --no-spliced-alignment \
         -q \
         -1 ${sra}_1.trim.fastq \
         -2 ${sra}_2.trim.fastq \
         -U ${sra}_1s.trim.fastq,${sra}_2s.trim.fastq \
         -S ${sra}_vs_${params.ref.prefix}.sam \
         -t \
         -p 1 \
         --dta-cufflinks
     else
       hisat2 \
         -x ${params.ref.prefix} \
         --no-spliced-alignment \
         -q \
         -U ${sra}_1.trim.fastq \
         -S ${sra}_vs_${params.ref.prefix}.sam \
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
  publishDir "$sra", mode: 'link'
  stageInMode "link"
  tag { sra }

  input:
    set val(sra), file("${sra}_vs_${params.ref.prefix}.sam") from INDEXED_SAMPLES

  output:
    set val(sra), file("${sra}_vs_${params.ref.prefix}.bam") into SORTED_FOR_INDEX

  script:
    """
    samtools sort -o ${sra}_vs_${params.ref.prefix}.bam -O bam ${sra}_vs_${params.ref.prefix}.sam
    """
}


/*
 * Indexes the BAM alignment file
 *
 * depends: samtools_index
 */
process samtools_index {
  module 'samtools'
  publishDir "$sra", mode: 'link'
  stageInMode "link"

  input:
    set val(sra), file("${sra}_vs_${params.ref.prefix}.bam") from SORTED_FOR_INDEX

  output:
    set val(sra), file("${sra}_vs_${params.ref.prefix}.bam"), file("${sra}_vs_${params.ref.prefix}.bam.bai") into BAM_INDEXED_FOR_STRINGTIE

  script:
    """
    samtools index ${sra}_vs_${params.ref.prefix}.bam
    """
}


/*
 * Generates expression-level transcript abundance
 *
 * depends: samtools_index
 */
process stringtie {
  module 'stringtie'
  publishDir "$sra", mode: 'link'
  stageInMode "link"

  input:
    // We don't really need the .bai file, but we want to ensure
    // this process runs after the samtools_index step so we
    // require it as an input file.
    set val(sra), file("${sra}_vs_${params.ref.prefix}.bam"), file("${sra}_vs_${params.ref.prefix}.bam.bai") from BAM_INDEXED_FOR_STRINGTIE

  output:
    set val(sra), file("${sra}_vs_${params.ref.prefix}.gtf") into STRINGTIE_GTF

  script:
    """
    stringtie -v -p 1 -e -G ${params.ref.path}${params.ref.prefix}.gtf -o ${sra}_vs_${params.ref.prefix}.gtf -l ${sra} ${sra}_vs_${params.ref.prefix}.bam
    """
}


/*
 * Generates the final FPKM file
 */
process fpkm {
  publishDir "$sra", mode: 'link'
  stageInMode "link"

  input:
    set val(sra), file("${sra}_vs_${params.ref.prefix}.gtf") from STRINGTIE_GTF

  output:
    file "${sra}_vs_${params.ref.prefix}.fpkm" into FPKMS

  script:
    """
    ${PWD}/scripts/gtf2fpkm.sh ${sra} ${params.ref.prefix}
    """
}
