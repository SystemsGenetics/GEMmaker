.. figure:: images/GEMmaker-logo-sm.png
   :alt: GEMmaker Logo

|DOI|

Software Prerequisites
----------------------


GEMmaker is a `Nextflow <https://www.nextflow.io/>`__ workflow for
large-scale gene expression sample processing, expression-level
quantification and Gene Expression Matrix (GEM) construction. Results
from GEMmaker are useful for differential gene expression (DGE) and gene
co-expression network (GCN) analyses. The GEMmaker workflow currently
supports Illumina RNA-seq datasets.

Before execution of GEMmaker you must have the necessary software. The
following list provides the set of tools and versions that have been
verified to work with GEMmaker. NOTE: newer versions of these tools are
assumed to also work, and older versions may work but have not been
tested:

-  `python3 <https://www.python.org>`__
-  `nextflow <https://www.nextflow.io/>`__ v0.32: Executes the workflow.
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
-  `samtools <http://www.htslib.org/>`__ v1.3.1: Used for indexing and
   sorting of BAM files created by Hisat2.
-  `stringTie <http://www.ccb.jhu.edu/software/stringtie/>`__ v1.3.4d:
   Performs gene expression quantification.
-  `MultiQC <http://multiqc.info/>`__ (optional) v1.5: Generate a full
   summary report for the entire workflow.

Additionally, GEMmaker requires that some of these packages be available
on your system as `Environment
Modules <http://modules.sourceforge.net/>`__ or
`Lmod <https://www.tacc.utexas.edu/research-development/tacc-projects/lmod>`__
modules and that these modules are named as follows:

-  fastQC
-  hisat2
-  sammtools
-  sratoolkit
-  stringtie
-  trimmomatic

Preparing the prerequisites may be challenging for some. Additionally,
some HPC systems may have the software available but not using the
module names listed above. There are a few options to simplify execution
of GEMmaker despite these problems.

Downloading GEMmaker
~~~~~~~~~~~~~~~~~~~~

After ensuring that all necessary software prerequisites are available,
clone GEMmaker into a working directory.

To clone the workflow into a directory:

.. code:: bash

    nextflow clone SystemsGenetics/GEMmaker target-dir

As with all Nextflow workflows, you can configure the behavior of the
workflow by creating a ``nextflow.config`` file. The GEMmaker workflow
provides an example file (``nextflow.config.example``) which you can
copy to get started.

.. code:: bash

    cp nextflow.config.example nextflow.config

Edit ``nextflow.config`` according to the instructions. At the
very least, you will need to modify the input parameters and the
execution profile. You may want to refer to the `Nextflow configuration
documentation <https://www.nextflow.io/docs/latest/config.html>`__ to
set proper profile settings for your environment. For example, to run
the workflow on an HPC system you will have to specify the "executor"
that corresponds to your system's scheduler (such as ``pbs``, ``slurm``,
etc), as well as any other properties specific to your system, such as
your job queue.


Local machine
~~~~~~~~~~~~~

You can execute GEMmaker out-of-the-box using the
`GEMmaker-docker <https://github.com/SystemsGenetics/GEMmaker-docker>`__
image. Refer to the Github repo for usage instructions. It contains all
of the necessary software needed to execute the workflow. No
installation of software dependencies is required and it will not
conflict with existing software. Note that execution of GEMmaker for a
large number of samples on a local machine is not recommended as it may
take a very long time to complete.

High-Performance Computing (HPC) cluster
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

To execute GEMmaker on an HPC cluster you must do **only one** of the
following: - Have your HPC admins install the necessary software using
the module names specified above. - Install the software into your own
space and create your own module files for the software. You can find
examples of module files in the ``files`` directory of the
`GEMmaker-docker <https://github.com/SystemsGenetics/GEMmaker-docker>`__
repository. - Edit ``main.nf`` and alter the module names to match those
of your HPC system.

.. |DOI| image:: https://zenodo.org/badge/114067776.svg
   :target: https://zenodo.org/badge/latestdoi/114067776
