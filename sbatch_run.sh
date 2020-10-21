#!/bin/sh
#SBATCH --partition=ficklin
#SBATCH --account=ficklin
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=7-00:00:00
#SBATCH --job-name=GalaApple
#SBATCH --output=GalaApple.out
#SBATCH --error=GalaApple.err
#SBATCH --mail-user=john.hadish@wsu.edu
#SBATCH --mail-type=ALL
#SBATCH --mem=256000 

module add java nextflow singularity/2.4.2
nextflow run main.nf -profile slurm,singularity -resume -with-report -with-timeline  -with-trace

