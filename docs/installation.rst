.. _installation:

Installation
------------

Dependencies
~~~~~~~~~~~~

GEMmaker requires a variety of common and bioinformatics software.  You must install the required software listed in the **Required Dependencies** section below.  Afterwards, there are two ways to ensure all other software dependencies are met:

1.  Use pre-built containerized versions of the software via `Docker <https://www.docker.com/>`__ or `Singularity <https://sylabs.io/>`__. **(Recommended)**
2.  Install all software manually, on your local machine or computational cluster and ensure that the software verions are compatible.


For easiest use, it is recommended to use pre-built containers. For users not familiar with containers or users without administrative control over the computational machine, it is recommended to use Singularity.  Example instructions provided in this documentation will assume Singularity is available.

Required Dependencies
*********************

At a minimum, GEMmaker requires the following:

- `java <https://www.java.com/en/>`__ v1.8.0 or later: Prerequisite of nextflow.
- `nextflow <https://www.nextflow.io/>`__ v18.10.1: Executes the workflow.
- `Singularity <https://sylabs.io/>`__ or `Docker <https://www.docker.com/>`__. **(Singularity Recommended)**

.. warning::
  Nextflow does not yet support Singularity 3.0 or later. You can use any version of Singularity between 2.4 and 2.6.

Container Support
*****************

All of the software tools needed to run GEMmaker have been pre-installed into containers by the GEMmaker development team. Therefore, you do not need to install them!  Using these containerized **images** can ensure that results from GEMmaker are always reproducible because the environment in which the software is executed will never change, even if the host computational computer is updated.  Please ensure one of these containerization services is installed.

  - `Singularity Community Edition (CE) <https://sylabs.io/>`__  **(recommended)**.
  - `Docker <https://www.docker.com/>`__.


.. warning::
  Nextflow does not yet support Singularity 3.0 or later. You can use any version of Singularity between 2.4 and 2.6.


Other Software Dependencies
***************************

The following software are required.  If you have opted to use containers with GEMmaker you **do not** need to install these. If you have opted for manual installation each of these tools must be available on the host machine.

.. note::
  If you opt for manual installation, other versions of these tools may work but have not been tested.


-  `python3 <https://www.python.org>`__ v3.7
-  `sratoolkit <https://www.ncbi.nlm.nih.gov/books/NBK158900/>`__ v2.10.0: Downloads SRA files from NCBI using the SRA Run IDs.
-  `fastQC <https://www.bioinformatics.babraham.ac.uk/projects/fastqc/>`__ v0.11.7: Generates read quality statistics for FASTQ files used by the workflow.
-  `trimmomatic <http://www.usadellab.org/cms/?page=trimmomatic>`__ v0.38: Removes low-quality bases from the ends of reads and removes adapter sequences.
-  `hisat2 <https://ccb.jhu.edu/software/hisat2/index.shtml>`__ v2.1.0: Aligns cleaned reads to the reference genome.
-  `salmon <https://combine-lab.github.io/salmon/>`__ v0.12.0: Performs quasi-alignment of reads and quantifies
-  `kallisto <https://pachterlab.github.io/kallisto/>`__ v0.45.0: Performs pseudo-alignment of reads and quantifies
-  `samtools <http://www.htslib.org/>`__ v1.3.1: Used for indexing and sorting of BAM files created by Hisat2.
-  `stringTie <http://www.ccb.jhu.edu/software/stringtie/>`__ v1.3.4d: Performs gene expression quantification.
-  `MultiQC <http://multiqc.info/>`__ v1.5: Generate a full summary report for the entire workflow.

Choice of Computing Environment
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Local Machine
*************

GEMmaker can be run on a local computer; in particular, it has been tested on `Ubuntu 18.04 <https://www.ubuntu.com/>`__.  GEMmaker has been optimized to allow for execution of large sample sizes without consuming large amounts of storage space.  However, the time required to execute is dependent on the amount of computing power the machine has. You should consider using a High-Performance Computing (HPC) system if you have a large number of samples and would like GEMmaker to execute faster.

High-Performance Computing (HPC) Setup
**************************************

On an HPC system it is recommended to use containerized dependencies as most HPC users do not have access to install software and HPC systems can change yielding results that may be hard to reproduce after time has passed and system setups have changed.  Additionally, most HPC setups do not allow users to run Docker, but rather provide Singularity instead. Using singularity is recommended on an HPC system rather than installing software dependencies manually. You will need to make sure that nextflow and Singularity are installed on your cluster (you may need the help of your HPC administrator).

If you do not have access to Singularity on your HPC, GEMmaker does support use of the `LMOD <https://lmod.readthedocs.io/en/latest/>`__ module system which is popular for managing software on HPCs. If you choose to install software manually, your HPC system admin will most likely need to install them for you.

Kubernetes
**********

GEMmaker can be used on a `Kubernetes <https://kubernetes.io/>`__ cluster with minimal effort. Consult the `kube-runner <https://github.com/SystemsGenetics/kube-runner>`__ project for instructions.

Installing GEMmaker
~~~~~~~~~~~~~~~~~~~

Once all of the necessary software dependencies have been installed, clone GEMmaker into your working directory and checkout the stable version of your choice:

.. code:: bash

  nextflow clone SystemsGenetics/GEMmaker
  cd GEMmaker
  checkout v1.1
