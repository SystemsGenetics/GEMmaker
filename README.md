[![DOI](https://zenodo.org/badge/114067776.svg)](https://zenodo.org/badge/latestdoi/114067776)
[![Build Status](https://travis-ci.org/SystemsGenetics/GEMmaker.svg?branch=master)](https://travis-ci.org/SystemsGenetics/GEMmaker)

![GEMmaker Logo](images/GEMmaker-logo-sm.png)

GEMmaker is a [Nextflow](https://www.nextflow.io/) workflow for large-scale gene expression sample processing, expression-level quantification and Gene Expression Matrix (GEM) construction. Results from GEMmaker are useful for differential gene expression (DGE) and gene co-expression network (GCN) analyses. This report is the MultiQC summary of the GEMmaker workflow results. The GEMmaker workflow currently supports Illumina RNA-seq datasets.

## Gemaker is:

### Easy to Use
![Ease of Use](images/ease_of_use.png)
1. No bioinformatics software installation required.
2. Runs on a stand-alone computer or High Performance Compute (HPC) cluster.
3. Simple configuration file setup.
4. Resulting data is ready for Differential Gene Expression (DEG) or Gene Co-Expression Network (GCN) analysis.
5. Full online documentation.

### Reproducible
![Reproducible](images/reproducible.png)
1. Software versions and computing environment are always the same every time GEMmaker is repeatedly run
2. Sharing input data and config files ensures anyone can reproduce exact results.

### Interoperable  
![Interoperable](images/interoperable.png)
1. Uses a variety of bioinformatics tools:  SRAToolkit, FastQC, Trimmomatic, Hisat2, Stringtie, Kallisto, Salmon, MultiQC
2. Integrates with iRODs for easy data movement
3. Easily retrieves samples from NCBI’s Sequence Read Archive (SRA)
4. Can combine local samples with those from SRA
5. Runs on any modern HPC scheduler.

### Findable
![Findable](images/findable_data.png)
1. Sample meta-data is retrieved from NCBI SRA
2. Controlled vocabularies are used to automatically remap SRA annotations
3. JSON-format meta-data files or each sample are provided
4. Meta-data files can be integration with databases and iRODs for querying and searching.

### Scalable
![Scalable](images/scalable.png)
1. Useful for small DGE projects of ~100 samples as well as large scale GCN projects with 1000s of samples.
2. Cleans unnecessary files as it goes.
3. Keeps storage requirements to a minimum.


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

## Acknowledgments
GEMmaker is a collaborative project of the [Ficklin](http://ficklinlab.cahnrs.wsu.edu/) and [Feltus](https://www.clemson.edu/science/departments/genetics-biochemistry/people/profiles/ffeltus) programs at [Washington State University](http://www.wsu.edu) and [Clemson University](http://www.clemson.edu) respectively with guidance from [RENCI](https://renci.org/).

GEMmaker is funded by the [NSF SciDAS](http://scidas.org/) project, [award #1659300](https://www.nsf.gov/awardsearch/showAward?AWD_ID=1659300)

!["WSU"](images/WSU.png) !["Clemson"](images/clemson.png)
!["RENCI"](images/renci.png)
!["NSF"](images/NSF.png)
