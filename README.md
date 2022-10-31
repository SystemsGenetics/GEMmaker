# ![GEMmaker](docs/images/GEMmaker-logo-sm.png)

**GEMmaker is a Nextflow workflow for large-scale gene expression sample processing, expression-level quantification and Gene Expression Matrix (GEM) construction. Results from GEMmaker are useful for differential gene expression (DGE) and gene co-expression network (GCN) analyses. The GEMmaker workflow currently supports Illumina RNA-seq datasets.**.

[![DOI](https://zenodo.org/badge/114067776.svg)](https://zenodo.org/badge/latestdoi/114067776)

[![GitHub Actions CI Status](https://github.com/systemsgenetics/gemmaker/workflows/nf-core%20CI/badge.svg)](https://github.com/systemsgenetics/gemmaker/actions?query=workflow%3A%22nf-core+CI%22)
[![GitHub Actions Linting Status](https://github.com/systemsgenetics/gemmaker/workflows/nf-core%20linting/badge.svg)](https://github.com/systemsgenetics/gemmaker/actions?query=workflow%3A%22nf-core+linting%22)

[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A521.04.0-23aa62.svg?labelColor=000000)](https://www.nextflow.io/)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)

## Introduction

GEMmaker (i.e. **systemsgenetics/gemmaker**) is a pipeline for quantification of Illumina RNA-seq data. Users can choose from Hisat2, Kallisto or Salmon. It can process locally stored FASTQ files or automatically retrieve them from NCBI's SRA.  The pipeline is built using [Nextflow](https://www.nextflow.io), a workflow tool to run tasks across multiple compute infrastructures in a very portable manner. It comes with docker containers making installation trivial and results highly reproducible.

## nf-core Compatibility

GEMmaker is an [nf-core](https://nf-co.re/) compatible workflow, however, GEMmaker is not an official nf-core workflow.  This is because nf-core offers the [nf-core/rnaseq](https://nf-co.re/rnaseq) workflow which is an excellent workflow for RNA-seq analysis that provides similar functionality to GEMmaker.  However, GEMmaker is different in that it can scale to thousands of samples without exceeding local storage resources by running samples in batches and removing intermediate files.  It can do the same for smaller sample sets on machines with less computational resources.  This ability to scale is a unique adaption that is currently not provided by Nextflow.   When Nextflow does provide support for batching and scaling, the [nf-core/rnaseq](https://nf-co.re/rnaseq) will be updated and GEMmaker will probably be retired in favor of the nf-core workflow. Until then, if you are limited by storage GEMmaker can help!
v

## How to Use

Please see the [GEMmaker documentation](https://gemmaker.readthedocs.io/en/latest/) for in-depth instructions for running GEMmaker.

## Credits

Please see the list of developers who have contributed to this repository.

Development of GEMmaker was funded by the U.S. National Science Foundation Award [#1659300](https://www.nsf.gov/awardsearch/showAward?AWD_ID=1659300&HistoricalAwards=false).

If you use GEMmaker in your research, please use this citation:

**Hadish, J. A., Biggs, T. D., Shealy, B. T., Bender, M. R., McKnight, C. B., Wytko, C., Smith, M. C., Feltus, F. A., Honaas, L., & Ficklin, S. P. (2022). GEMmaker: process massive RNA-seq datasets on heterogeneous computational infrastructure. BMC Bioinformatics, 23(1), 1–11.**

## Contributions and Support

If you would like to contribute to this pipeline, please see the [contributing guidelines](.github/CONTRIBUTING.md).

## Quick Start

Please follow the instructions in the ['Online Documentation'](https://gemmaker.readthedocs.io/en/latest/)
