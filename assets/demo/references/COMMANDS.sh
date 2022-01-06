#! /bin/bash

#This is a command to make the hisat files
module load hisat2
hisat2-build -f CORG.fna CORG | tee > hisat2-build.log

# The testing of salmon and KALLISTO was done in this directory. I need to find a better conversion toool than TOPHAT for .gtf + .fna conversion to transcripts.
# Also, I need to find a way to do gene abundance with it (currently it is transcript based, which is fine for CORG, but not for Arabidopsis for example)
#Tophat command to get fasta transcript file from the fna and the gtf file:
module load tophat
gtf_to_fasta CORG.gtf CORG.fna CORG.transcripts
# tophat assigns numbers as the first, so we have to get the gene name
# Produces:
#  >0 Alpha.1 Chr1. 1-100
#  AGTTTGGCCCGTGAAAGAAAGAAAAACAAAACTTATTTATTGAAAAATGACACGTGCAGA
#  ATCATTTACTAAACCAGACAATAACCGGTTAACGGTTTAA
#
# But we want:
#  >Alpha.1 Chr1. 1-100
#  AGTTTGGCCCGTGAAAGAAAGAAAAACAAAACTTATTTATTGAAAAATGACACGTGCAGA
#  ATCATTTACTAAACCAGACAATAACCGGTTAACGGTTTAA
# So we use the following sed to "clean" the top.hat output
cat CORG.transcripts | sed -r 's/(>[0-9]*\s)(.*)/>\2/g' > CORG.transcripts.clean
# Of course, all od this can be avoided if you just download a version of all the genes you want to align to...
# Here is Kallisto's Perl script to do the same thing I did with a sed. I could not get the perl to work: https://github.com/griffithlab/rnaseq_tutorial/wiki/Kallisto

#These are the commands for making the Kalisto Index:
module load kallisto
kallisto index -i CORG.transcripts.Kalisto.indexed CORG.transcripts.clean

#This is the commands for making the Salmon Index:
module load salmon
salmon index -t CORG.transcripts.clean -i CORG.transcripts.Salmon.indexed

#These are the commands for making the Star Index:
mkdir CORG.genome.Star.indexed

singularity exec -B ${PWD} https://depot.galaxyproject.org/singularity/star:2.7.9a--h9ee0642_0 STAR \
--runThreadN 4 \
--runMode genomeGenerate \
--genomeDir CORG.genome.Star.indexed \
--genomeFastaFiles CORG.fna \
--sjdbGTFfile CORG.transcripts.gtf
