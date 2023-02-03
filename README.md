![GEMmaker](docs/images/GEMmaker-logo-sm.png)

[![DOI](https://zenodo.org/badge/114067776.svg)](https://zenodo.org/badge/latestdoi/114067776)

[![GitHub Actions CI Status](https://github.com/systemsgenetics/gemmaker/workflows/nf-core%20CI/badge.svg)](https://github.com/systemsgenetics/gemmaker/actions?query=workflow%3A%22nf-core+CI%22)
[![GitHub Actions Linting Status](https://github.com/systemsgenetics/gemmaker/workflows/nf-core%20linting/badge.svg)](https://github.com/systemsgenetics/gemmaker/actions?query=workflow%3A%22nf-core+linting%22)

[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A521.04.0-23aa62.svg?labelColor=000000)](https://www.nextflow.io/)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)

## Introduction

**systemsgenetics/gemmaker** is a bioinformatics best-practice analysis pipeline for A workflow for construction of Gene Expression count Matrices (GEMs). Useful for Differential Gene Expression (DGE) analysis and Gene Co-E
xpression Network (GCN) construction
.

The pipeline is built using [Nextflow](https://www.nextflow.io), a workflow tool to run tasks across multiple compute infrastructures in a very portable manner. It uses Docker/Singularity containers making installation trivial and results highly reproducible. The [Nextflow DSL2](https://www.nextflow.io/docs/latest/dsl2.html) implementation of this pipeline uses one container per process which makes it much easier to maintain and update software dependencies. Where possible, these processes have been submitted to and installed from [nf-core/modules](https://github.com/nf-core/modules) in order to make them available to all nf-core pipelines, and to everyone within the Nextflow community!


On release, automated continuous integration tests run the pipeline on a full-sized dataset on the AWS cloud infrastructure. This ensures that the pipeline runs on AWS, has sensible resource allocation defaults set to run on real-world datasets, and permits the persistent storage of results to benchmark between pipeline releases and other analysis sources.


## Quick Start

Please see the [GEMmaker documentation](https://gemmaker.readthedocs.io/en/latest/) for in-depth instructions for running GEMmaker.

## Credits

GEMmaker was originally written by John Hadish, Tyler Biggs, Ben Shealy, Connor Wytko, Sai Prudhvi Oruganti, F. Alex Feltus, & Stephen Ficklin.

Development of GEMmaker was funded by the U.S. National Science Foundation Award [#1659300](https://www.nsf.gov/awardsearch/showAward?AWD_ID=1659300&HistoricalAwards=false).

## Contributions and Support

If you would like to contribute to this pipeline, please see the [contributing guidelines](.github/CONTRIBUTING.md).

## Citations

An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.

This pipeline uses code and infrastructure developed and maintained by the [nf-core](https://nf-co.re) community, reused here under the [MIT license](https://github.com/nf-core/tools/blob/master/LICENSE).

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).
