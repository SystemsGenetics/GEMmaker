# ![GEMmaker](docs/images/GEMmaker-logo-sm.png)

**GEMmaker is a Nextflow workflow for large-scale gene expression sample processing, expression-level quantification and Gene Expression Matrix (GEM) construction. Results from GEMmaker are useful for differential gene expression (DGE) and gene co-expression network (GCN) analyses. The GEMmaker workflow currently supports Illumina RNA-seq datasets.**.

[![DOI](https://zenodo.org/badge/114067776.svg)](https://zenodo.org/badge/latestdoi/114067776)
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

### GEMmaker
If you use GEMmaker for your analysis, please cite it using the following:

> John Hadish, Tyler Biggs, Ben Shealy, Connor Wytko, Sai Prudhvi Oruganti, F. Alex Feltus, & Stephen Ficklin. (2020, January 22). SystemsGenetics/GEMmaker: Release v1.1 (Version v1.1). Zenodo. http://doi.org/10.5281/zenodo.3620945
You can cite the `nf-core` publication as follows:

### nf-core
To cite the nf-core framework for community-curated bioinformatics pipelines:
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.5281/zenodo.3620945](http://doi.org/10.5281/zenodo.3620945).

### Tools used by GEMmaker
In addition, references or URLs for tools used in this pipeline are as follows:

#### Nextflow
> Di Tommaso, P., Chatzou, M., Floden, E. W., Barja, P. P., Palumbo, E., & Notredame, C. (2017). Nextflow enables reproducible computational workflows. Nature Biotechnology, 35(4), 316–319. https://doi.org/10.1038/nbt.3820

#### SRAtoolkit
https://github.com/ncbi/sra-tools

#### Aspera
https://www.ibm.com/products/aspera

#### FastQC
https://www.bioinformatics.babraham.ac.uk/projects/fastqc/

#### Trimmomatic
> Bolger, A. M., Lohse, M., & Usadel, B. (2014). Trimmomatic: A flexible trimmer for Illumina sequence data. Bioinformatics, 30(15), 2114–2120. https://doi.org/10.1093/bioinformatics/btu170

#### Hisat2
> Kim, D., Paggi, J. M., Park, C., Bennett, C., & Salzberg, S. L. (2019). Graph-based genome alignment and genotyping with HISAT2 and HISAT-genotype. Nature Biotechnology, 37(8), 907–915. https://doi.org/10.1038/s41587-019-0201-4

#### Kallisto
> Bray, N. L., Pimentel, H., Melsted, P., & Pachter, L. (2016). Near-optimal probabilistic RNA-seq quantification. Nature Biotechnology. https://doi.org/10.1038/nbt.3519

#### Salmon
> Patro, R., Duggal, G., Love, M. I., Irizarry, R. A., & Kingsford, C. (2017). Salmon provides fast and bias-aware quantification of transcript expression. Nature Methods, 14(4), 417–419. https://doi.org/10.1038/nmeth.4197

#### SAMtools
> Li, H., Handsaker, B., Wysoker, A., Fennell, T., Ruan, J., Homer, N., Marth, G., Abecasis, G., & Durbin, R. (2009). The Sequence Alignment/Map format and SAMtools. Bioinformatics, 25(16), 2078–2079. https://doi.org/10.1093/bioinformatics/btp352

#### StringTie
> Pertea, M., Kim, D., Pertea, G. M., Leek, J. T., & Salzberg, S. L. (2016). Transcript-level expression analysis of RNA-seq experiments with HISAT, StringTie and Ballgown. Nature Protocols, 11(9), 1650–1667. https://doi.org/10.1038/nprot.2016.095

#### MultiQC
> Pertea, M., Kim, D., Pertea, G. M., Leek, J. T., & Salzberg, S. L. (2016). Transcript-level expression analysis of RNA-seq experiments with HISAT, StringTie and Ballgown. Nature Protocols, 11(9), 1650–1667. https://doi.org/10.1038/nprot.2016.095
