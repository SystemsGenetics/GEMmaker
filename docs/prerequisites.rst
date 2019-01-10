.. figure:: images/GEMmaker-logo-sm.png
   :alt: GEMmaker Logo

|DOI|

.. _software_prerequisites:

Software Prerequisites
----------------------


GEMmaker is a `Nextflow <https://www.nextflow.io/>`__ workflow for large-scale
gene expression sample processing, expression-level quantification and Gene
Expression Matrix (GEM) construction. Results from GEMmaker are useful for
differential gene expression (DGE) and gene co-expression network (GCN)
analyses. The GEMmaker workflow currently supports Illumina RNA-seq datasets.
GEMmaker is intended to be run on a linux based operating system, the defacto
operating system of bioinformatics.

GEMmaker is compatible with `Docker <https://www.docker.com/>`__ and
`Singularity <https://www.sylabs.io/docs/>`__. **Running with either Docker or
Singularity is the recommended way to run GEMmaker**, as GEMmaker will
automatically pull each of the images that it needs to run.

.. note::

  Nextflow (The workflow language GEMmaker uses) does not support Singualrity
  version 3.0 and higher as of yet. You must use a version of singularity that
  is higher than 2.3.x, but less than 3.0.x

The following list provides the set of software and versions that have been
verified to work with GEMmaker. They are split into 2 sections. **Section 1**
are prerequisites of GEMmaker, and must be installed manually by the user
previous to running GEMmaker. **Section 2** are software contained in Docker
images that GEMmaker will automatically download. Most users of GEMmaker do not
need to worry about programs in section 2 unless a different version of the
software is required. NOTE: newer versions of these tools are assumed to also
work, and older versions may work but have not been tested:

**Section 1**: Required Prerequisites:

-  `java <https://www.java.com/en/>`__ v8 or later: prerequisite of nextflow
-  `nextflow <https://www.nextflow.io/>`__ v18.10.1: Executes the workflow. Must
   be installed on system by user.
-  either `Docker <https://www.docker.com/>`__
   v18.09.0 or `Singularity <https://www.sylabs.io/docs/>`__ v2.5 - 2.6 NOT
   v3.0!: Automatically downloads and runs the docker images of the following
   programs:

**Section 2**: Programs contained in Docker Images that will be automatically
downloaded by GEMmaker using Docker or Singularity:

-  `python3 <https://www.python.org>`__
-  `sratoolkit <https://www.ncbi.nlm.nih.gov/books/NBK158900/>`__
   v2.8.0: Downloads SRA files from NCBI using the SRA Run IDs.
-  `fastQC <https://www.bioinformatics.babraham.ac.uk/projects/fastqc/>`__
   v0.11.7: Generates read quality statistics for FASTQ files used by
   the workflow.
-  `trimmomatic <http://www.usadellab.org/cms/?page=trimmomatic>`__
   v0.38: Removes low-quality bases from the ends of reads and removes
   adapter sequences.
-  `hisat2 <https://ccb.jhu.edu/software/hisat2/index.shtml>`__ v2.1.0:
   Aligns cleaned reads to the reference genome.
-  `salmon <https://combine-lab.github.io/salmon/>`__ v0.12.0:
   Performs quasi alignment of reads and quantifies
-  `kallisto <https://pachterlab.github.io/kallisto/>`__ v 0.45.0
   Performs pseudo alignment of reads and quantifies
-  `samtools <http://www.htslib.org/>`__ v1.3.1: Used for indexing and
   sorting of BAM files created by Hisat2.
-  `stringTie <http://www.ccb.jhu.edu/software/stringtie/>`__ v1.3.4d:
   Performs gene expression quantification.
-  `MultiQC <http://multiqc.info/>`__ (optional) v1.5: Generate a full
   summary report for the entire workflow.

If a user requires a different version of software contained in **Section 2**,
please reference **different software versions**

Local machine
~~~~~~~~~~~~~

GEMmaker can be ran on a local desktop machine, assuming that the machine is
running linux, (we like `ubuntu 18.04 <https://www.ubuntu.com/>`__) has a
reasonable amount of  memory and has docker or singularity installed.


High-Performance Computing (HPC) cluster
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Most HPC clusters do not allow users to run docker, but singularity should be
allowed. On your own, or with the help of your HPC administrator, install
singularity and nextflow if they are not already installed. After Nextflow and
singuarity are installed, clone GEMmaker into your desired directory as shown
. GEMmaker will then install all of the dependencies the first time it is
run.

Downloading GEMmaker
~~~~~~~~~~~~~~~~~~~~

After ensuring that all necessary software prerequisites from **Section 1** are
installed, clone GEMmaker into your working directory.

To clone the workflow into a directory:

.. code:: bash

    nextflow clone SystemsGenetics/GEMmaker /target-dir/

Note that you should change ``target-dir`` to your target directory.

As with all Nextflow workflows, you can configure the behavior of the workflow
by creating a ``nextflow.config`` file. The GEMmaker workflow provides an
example file (``nextflow.config.example``) which you can copy to get started:

.. code:: bash

    cp nextflow.config.example nextflow.config

After Installation
~~~~~~~~~~~~~~~~~~

**To run the example data, refer to** :ref:`running_the_examples`

**To get started running your own data, refer to** :ref:`running_your_data`

**To learn about all the parameters in the** ``nextflow.config`` **refer to the instructions
at:** :ref:`nextflow_config_instructions`.



.. |DOI| image:: https://zenodo.org/badge/114067776.svg
   :target: https://zenodo.org/badge/latestdoi/114067776
