.. figure:: images/GEMmaker-logo-sm.png
   :alt: GEMmaker Logo

|DOI|

.. _software_prerequisites:

Software Prerequisites
----------------------

GEMmaker is compatible with `Docker <https://www.docker.com/>`__ and
`Singularity <https://www.sylabs.io/docs/>`__. **Running with either Docker or
Singularity is the recommended way to run GEMmaker**, as GEMmaker will
automatically pull each of the images that it needs to run.

.. warning::

  Nextflow (the workflow language that GEMmaker uses) does not support Singularity 3.0
  or later as of yet. You must use Singularity version 2.4 or later, but earlier than 3.0.

The following list provides the set of software and versions that have been
verified to work with GEMmaker. They are split into two sections. **Section 1**
lists the prerequisites of GEMmaker, which must be installed manually by the user
before running GEMmaker. **Section 2** lists the software contained in Docker
images that GEMmaker will automatically download. You most likely do not
need to worry about programs in Section 2 unless you would like to run GEMmaker
without using Docker or Singularity.

.. note::
  Other versions of these tools may work but have not been tested:

**Section 1: Required Prerequisites**

-  `java <https://www.java.com/en/>`__ v1.80 or later: Prerequisite of nextflow.
-  `nextflow <https://www.nextflow.io/>`__ v18.10.1: Executes the workflow.
-  `Docker <https://www.docker.com/>`__ v18.09.0: Automatically downloads and runs the containerized dependencies.
-  `Singularity <https://www.sylabs.io/docs/>`__ v2.4 - v2.6: Alternative to Docker on systems where Docker is not available.

**Section 2: Containerized Prerequisites**

-  `python3 <https://www.python.org>`__
-  `sratoolkit <https://www.ncbi.nlm.nih.gov/books/NBK158900/>`__ v2.9.2: Downloads SRA files from NCBI using the SRA Run IDs.
-  `fastQC <https://www.bioinformatics.babraham.ac.uk/projects/fastqc/>`__ v0.11.7: Generates read quality statistics for FASTQ files used by the workflow.
-  `trimmomatic <http://www.usadellab.org/cms/?page=trimmomatic>`__ v0.38: Removes low-quality bases from the ends of reads and removes adapter sequences.
-  `hisat2 <https://ccb.jhu.edu/software/hisat2/index.shtml>`__ v2.1.0: Aligns cleaned reads to the reference genome.
-  `salmon <https://combine-lab.github.io/salmon/>`__ v0.12.0: Performs quasi-alignment of reads and quantifies
-  `kallisto <https://pachterlab.github.io/kallisto/>`__ v0.45.0: Performs pseudo-alignment of reads and quantifies
-  `samtools <http://www.htslib.org/>`__ v1.3.1: Used for indexing and sorting of BAM files created by Hisat2.
-  `stringTie <http://www.ccb.jhu.edu/software/stringtie/>`__ v1.3.4d: Performs gene expression quantification.
-  `MultiQC <http://multiqc.info/>`__ (optional) v1.5: Generate a full summary report for the entire workflow.

Running on a stand-alone machine
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

GEMmaker can be ran on a local stand-alone desktop machine, assuming that the machine is
running linux, (specifically it has been tested on `ubuntu 18.04 <https://www.ubuntu.com/>`__) has a
reasonable amount of RAM and has docker or singularity installed.

.. note::

  Generally speaking, local machines will not scale well to experiments larger than 10s of samples due to the amount of storage and CPU cores needed for large experiments.

Running on a High-Performance Computing (HPC) cluster
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Most HPC clusters do not allow users to run Docker, but Singularity is often allowed. On your own, if you have permissions, or with the help of your HPC administrator, install Singularity and nextflow if they are not already installed. After Nextflow and Singularity are installed, clone GEMmaker into your desired directory as shown. GEMmaker will then install all of the dependencies the first time it is run.

If your HPC environment does not support Singularity, you can still run GEMmaker by installing the containerized dependencies as Environment Modules. You may need to work with your HPC administrator to install the necessary software and make it available via the module system.

Downloading GEMmaker
~~~~~~~~~~~~~~~~~~~~

After ensuring that all necessary software prerequisites from **Section 1** are
installed, clone GEMmaker into your working directory.

To clone the workflow into a directory:

.. code:: bash

    nextflow clone SystemsGenetics/GEMmaker target-dir/

Note that you should change ``target-dir`` to your target directory.

As with all Nextflow workflows, you should configure the behavior of the workflow by creating a ``nextflow.config`` file. The GEMmaker workflow provides an example file (``nextflow.config.example``) which you can copy to get started:

.. code:: bash

    cp nextflow.config.example nextflow.config

The example workflow is set up to execute example data that is provided in the GEMmaker repo.

After Installation
~~~~~~~~~~~~~~~~~~

**To run the example data, refer to** :ref:`running_the_examples`

**To get started running your own data, refer to** :ref:`running_your_data`

**To learn about all the parameters in the** ``nextflow.config`` **refer to the instructions
at:** :ref:`nextflow_config_instructions`.



.. |DOI| image:: https://zenodo.org/badge/114067776.svg
   :target: https://zenodo.org/badge/latestdoi/114067776
