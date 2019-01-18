.. _installation:

Installation
------------

Dependencies
~~~~~~~~~~~~

GEMmaker depends primarily on `Nextflow <https://www.nextflow.io/>`__ for workflow execution and `Docker <https://www.docker.com/>`__ / `Singularity <https://www.sylabs.io/docs/>`__ for running other software dependencies in containers. **It is highly recommended that you use GEMmaker with Docker or Singularity**, as GEMmaker will automatically pull any other dependencies that it needs to run.

.. warning::
  Nextflow does not yet support Singularity 3.0 or later. You can use any version of Singularity between 2.4 and 2.6.

The following sections lists the software and versions that have been tested with GEMmaker. **Required Dependencies** are the dependencies which must be installed by the user before running GEMmaker. **Containerized Dependencies** are the dependencies which are contained in Docker images, which GEMmaker will automatically download. You most likely do not need to worry about the containerized dependencies unless you would like to run GEMmaker without using Docker or Singularity.

.. note::
  Other versions of these tools may work but have not been tested:

**Required Dependencies**

-  `java <https://www.java.com/en/>`__ v1.8.0 or later: Prerequisite of nextflow.
-  `nextflow <https://www.nextflow.io/>`__ v18.10.1: Executes the workflow.
-  `Docker <https://www.docker.com/>`__ v18.09.0: Automatically downloads and runs the containerized dependencies.
-  `Singularity <https://www.sylabs.io/docs/>`__ v2.4 - v2.6: Alternative to Docker on systems where Docker is not available.

**Containerized Dependencies**

-  `python3 <https://www.python.org>`__ v3.5.1
-  `sratoolkit <https://www.ncbi.nlm.nih.gov/books/NBK158900/>`__ v2.9.2: Downloads SRA files from NCBI using the SRA Run IDs.
-  `fastQC <https://www.bioinformatics.babraham.ac.uk/projects/fastqc/>`__ v0.11.7: Generates read quality statistics for FASTQ files used by the workflow.
-  `trimmomatic <http://www.usadellab.org/cms/?page=trimmomatic>`__ v0.38: Removes low-quality bases from the ends of reads and removes adapter sequences.
-  `hisat2 <https://ccb.jhu.edu/software/hisat2/index.shtml>`__ v2.1.0: Aligns cleaned reads to the reference genome.
-  `salmon <https://combine-lab.github.io/salmon/>`__ v0.12.0: Performs quasi-alignment of reads and quantifies
-  `kallisto <https://pachterlab.github.io/kallisto/>`__ v0.45.0: Performs pseudo-alignment of reads and quantifies
-  `samtools <http://www.htslib.org/>`__ v1.3.1: Used for indexing and sorting of BAM files created by Hisat2.
-  `stringTie <http://www.ccb.jhu.edu/software/stringtie/>`__ v1.3.4d: Performs gene expression quantification.
-  `MultiQC <http://multiqc.info/>`__ v1.5: Generate a full summary report for the entire workflow.

Environments
~~~~~~~~~~~~

Local machine
=============

GEMmaker can be run on a local computer; in particular, it has been tested on `Ubuntu 16.04 <https://www.ubuntu.com/>`__. However, local machines in general will not scale well to experiments larger than 10s of samples due to the high compute and storage requirements of large experiments.

High-Performance Computing (HPC) cluster
========================================

Most HPC clusters do not allow users to run Docker, but provide Singularity instead. You will need to make sure that nextflow and Singularity are installed on your cluster (you may need the help of your HPC administrator). If your HPC environment does not support Singularity, you will need to install the containerized dependencies as Environment Modules. Again, you may need to work with your HPC administrator to install the necessary software and make it available via the module system.

Installing GEMmaker
~~~~~~~~~~~~~~~~~~~

Once all of the necessary software dependencies have been installed, clone GEMmaker into your working directory:

.. code:: bash

  nextflow clone SystemsGenetics/GEMmaker
