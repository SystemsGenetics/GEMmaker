# ![GEMmaker](docs/images/GEMmaker-logo-sm.png)

**GEMmaker is a Nextflow workflow for large-scale gene expression sample processing, expression-level quantification and Gene Expression Matrix (GEM) construction. Results from GEMmaker are useful for differential gene expression (DGE) and gene co-expression network (GCN) analyses. The GEMmaker workflow currently supports Illumina RNA-seq datasets.**.

[![GitHub Actions CI Status](https://github.com/systemsgenetics/gemmaker/workflows/nf-core%20CI/badge.svg)](https://github.com/systemsgenetics/gemmaker/actions)
[![GitHub Actions Linting Status](https://github.com/systemsgenetics/gemmaker/workflows/nf-core%20linting/badge.svg)](https://github.com/systemsgenetics/gemmaker/actions)
[![Nextflow](https://img.shields.io/badge/nextflow-%E2%89%A520.04.0-brightgreen.svg)](https://www.nextflow.io/)

[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg)](https://bioconda.github.io/)
[![Docker](https://img.shields.io/docker/automated/systemsgenetics/gemmaker.svg)](https://hub.docker.com/r/systemsgenetics/gemmaker)

## Introduction

GEMmaker (i.e. **systemsgenetics/gemmaker**) is a pipeline for quantification of Illumina RNA-seq data. Users can choose from Hisat2, Kallisto or Salmon. It can process locally stored FASTQ files or automatically retrieve them from NCBI's SRA.  The pipeline is built using [Nextflow](https://www.nextflow.io), a workflow tool to run tasks across multiple compute infrastructures in a very portable manner. It comes with docker containers making installation trivial and results highly reproducible.

## nf-core Compatibility
GEMmaker is an [nf-core](https://nf-co.re/) compatible workflow, however, GEMmaker is not an official nf-core workflow.  This is because nf-core offers the [nf-core/rnaseq](https://nf-co.re/rnaseq) workflow which is an excellent workflow for RNA-seq analysis that provides similar functionality to GEMmaker.  However, GEMmaker is different in that it can scale to thousands of samples without exceeding local storage resources by running samples in batches and removing intermediate files.  It can do the same for smaller sample sets on machines with less computational resources.  This ability to scale is a unique adaption that is currently not provided by Nextflow.   When Nextflow does provide support for batching and scaling, the [nf-core/rnaseq](https://nf-co.re/rnaseq) will be updated and GEMmaker will probably be retired in favor of the nf-core workflow. Until then, if you are limited by storage GEMmaker can help!
v

## How to Use

Please see the [GEMmaker documentation](https://gemmaker.readthedocs.io/en/latest/) for in-depth instructions for running GEMmaker.


## Credits

GEMmaker was originally written by John Hadish, Tyler Biggs, Ben Shealy, Connor Wytko, Sai Prudhvi Oruganti, F. Alex Feltus, & Stephen Ficklin.


## Contributions and Support

If you would like to contribute to this pipeline, please see the [contributing guidelines](.github/CONTRIBUTING.md).

## Citations

If you use GEMmaker for your analysis, please cite it using the following doi: [10.5281/zenodo.XXXXXX](https://doi.org/10.5281/zenodo.XXXXXX)

You can cite the `nf-core` publication as follows:

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.5281/zenodo.3620945](http://doi.org/10.5281/zenodo.3620945).

In addition, references of tools and data used in this pipeline are as follows:

<!-- TODO nf-core: Add bibliography of tools and data used in your pipeline -->
