[![DOI](https://zenodo.org/badge/114067776.svg)](https://zenodo.org/badge/latestdoi/114067776)


![GEMmaker Logo](images/GEMmaker-logo-sm.png)


GEMmaker is a [Nextflow](https://www.nextflow.io/) workflow for large-scale gene expression sample processing, expression-level quantification and Gene Expression Matrix (GEM) construction. Results from GEMmaker are useful for differential gene expression (DGE) and gene co-expression network (GCN) analyses. This report is the MultiQC summary of the GEMmaker workflow results. The GEMmaker workflow currently supports Illumina RNA-seq datasets.


## Tools

GEMmaker uses the following tools:

- [python3](https://www.python.org)
- [nextflow](https://www.nextflow.io/) v0.32:  Executes the workflow.
- [sratoolkit](https://www.ncbi.nlm.nih.gov/books/NBK158900/) v2.8.0:  Downloads SRA files from NCBI using the SRA Run IDs.
- [fastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) v0.11.7:  Generates read quality statistics for FASTQ files used by the workflow.
- [trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic) v0.38:  Removes low-quality bases from the ends of reads and removes adapter sequences.
- [hisat2](https://ccb.jhu.edu/software/hisat2/index.shtml) v2.1.0:  Aligns cleaned reads to the reference genome.
- [samtools](http://www.htslib.org/) v1.3.1:  Used for indexing and sorting of BAM files created by Hisat2.
- [stringTie](http://www.ccb.jhu.edu/software/stringtie/) v1.3.4d:  Performs gene expression quantification.
- [MultiQC](http://multiqc.info/) (optional) v1.5:  Generate a full summary report for the entire workflow.

## Installation 
See the online GEMmaker [Documentation](https://gemmaker.readthedocs.io/en/latest/)

## Execution 
See the online GEMmaker [Documentation](https://gemmaker.readthedocs.io/en/latest/)
