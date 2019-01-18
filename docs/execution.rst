.. _execution:

Workflow Execution
------------------

Acquiring RNA-seq Samples
~~~~~~~~~~~~~~~~~~~~~~~~~

GEMmaker supports processing of sample files from your local computer and from the `NCBI SRA repository <https://www.ncbi.nlm.nih.gov/sra>`__.

-  For local samples, identify a `glob pattern <https://en.wikipedia.org/wiki/Glob_(programming)>`__ that finds these files.
-  For remote samples, identify the NCBI SRA run IDs of the samples you want to process. The run IDs typically start with an SRR, ERR, or DRR prefix. These run IDs must be placed, one per line, in a file and the filename should be set in ``remote_list_path`` of ``nextflow.config``.

Acquiring Reference Genome Files
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Download the genome annotation/reference files. These will differ depending on your organism that you are considering. The necessary files will also be different  depending on whether you use Hisat2, Salmon or Kallisto. The following instructions will describe how to make the proper files for each of these programs.

Hisat2
======

You must aquire two files from the Internet:

-  A FASTA file containing the full genomic sequence (either pseudomolecules or scaffolds). Note, if your genome file is extremely large with hundreds of thousands of contigs/scaffolds, you may want to reduce the size of the FASTA file to contain only those contigs/scaffolds with predicted annotated genes. This will be converted into the hisat2 index files using the commands below.
-  A `GTF <https://uswest.ensembl.org/info/website/upload/gff.html>`__ file containing the gene models.

.. note::
  Sometimes a genome assembly does not provide a GTF file, but rather provides a `GFF3 <https://uswest.ensembl.org/info/website/upload/gff.html>`__ file. You can convert the GFF file to a GTF file using the ``gffread`` tool from `cufflinks <http://cole-trapnell-lab.github.io/cufflinks/file_formats/>`__, which you may have to download separately. Here is an example command-line to convert a GFF3 to GTF:

  .. code:: bash

    gffread <gff_file> -T -o <gtf_file>

  The arguments ``<gff_file>`` and ``<gtf_file>`` should be substituted for the names of your GFF3 and desired GTF file respectively.

Additional considerations:

-  You must have hisat2 index files of your genome sequence. These are constructed by using the ``hisat2-build`` command.

  .. code:: bash

    hisat2-build -f YourFastaFile.fna YourPrefix | tee > hisat2-build.log

-  The GTF file and the hisat2 index files must have the same prefix and this prefix must be specified by ``reference_prefix`` in ``nextflow.config``.

  .. code:: bash

    CORG.1.ht2
    CORG.2.ht2
    CORG.3.ht2
    CORG.4.ht2
    CORG.5.ht2
    CORG.6.ht2
    CORG.7.ht2
    CORG.8.ht2
    CORG.fna
    CORG.gtf

-  All of the genome annotation files must be in a directory and this directory must be specified by ``reference_path`` in ``nextflow.config``.

Salmon and Kallisto
===================

The process for these two aligners is rather similar.

You must aquire this file from the Internet:

- A FASTA file containing the individual transcripts of the organism of interest. This should be formatted so that the name of each transcript is before the sequence of its corresponding transcipt. This is different from the FASTA file used for HISAT2.

This file needs to be converted into the proper index file for either Salmon or Kallistou using the following commands (assume the name of the transcript FASTA file is ``YourFasta.transcripts.fna``):

For Salmon:

.. code:: bash

  salmon index -t YourFasta.transcripts.fna -i YourFasta.transcripts.Salmon.indexed

For Kallisto:

.. code:: bash

  kallisto index -i YourFasta.transcripts.Kalisto.indexed YourFasta.transcripts.fna

.. note::
  If you are running GEMmaker with Docker images, you will have to run these commands from within the corresponding Docker image:

  .. code:: bash

    # with docker
    docker run --rm -it systemsgenetics/hisat:2.1.0 bash

    # with singularity
    singularity shell work-singularity/systemsgenetics-hisat2-2.1.0.img

Executing the Workflow
~~~~~~~~~~~~~~~~~~~~~~

To execute the workflow on a local machine:

.. code:: bash

  nextflow run main.nf -profile standard

To resume a workflow in the event of a failure:

.. code:: bash

  nextflow run main.nf -profile standard -resume

To execute the workflow and generate trace, timeline and execution reports:

.. code:: bash

  nextflow run main.nf -profile standard -with-report -with-timeline -with-trace

To execute the workflow on an HPC system you must edit ``nextflow.config`` and add an appropriate profile for your system. Refer to the `Nextflow documentation <https://www.nextflow.io/docs/latest/config.html#config-profiles>`__. You can then use any of the above commands by changing the ``-profile`` argument to use your profile.

Performance Considerations
==========================

For large experiments on an HPC system, it is important to make sure that you are effectively utilizing the resources of the system. There are a number of parameters in ``nextflow.config`` which can be used to increase performance based on the capabilities of your system:

- ``params.execution.threads``: All processes which support multithreading (such as trimmomatic) will use this number of threads. This setting should be determined by the number of cores per node on your system; for example, if your system has nodes with 16 cores per node then you could set the number of threads to 16 to make full use of those nodes.

- ``params.execution.queue_size``: Nextflow will only run up to 100 processes at a time by default, but you may be able to increase this value based on the queue limits of your system.

Generating a Summary Report
===========================

The `MultiQC <http://multiqc.info>`__ tool will automatically generate a report on how each process ran.

Generating a Gene Expression Matrix (GEM)
=========================================

After GEMmaker completes, the resulting GEMs will be output to ``output/GEMs/`` by default. This directory contains the final gene-expression matrices in raw, TPM and FPKM form, depending on which output formats are enabled in ``nextflow.config``.

Using GEMs in Other Workflows
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

DGE Analysis
============

The raw GEM can be used for DGE analysis in edgeR and other DGE software.

Network Analysis
================

Any GEM can be used to construct a gene-coexpression network (GCN). `KINC <https://github.com/SystemsGenetics/KINC>`__ is a high-performance application that can construct networks using Pearson or Spearman for pairwise correlation, as well as Gassian mixture models (GMMs) for pairwise clustering. KINC is a Qt/`ACE <https://github.com/SystemsGenetics/ACE>`__ application that is capable of running on CPUs and GPUs, which means that it can scale to larger workloads.
