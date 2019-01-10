.. _running_your_data:

Running your Data
-----------------


Step 1) RNA-seq sample location
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


GEMmaker supports processing of sample files that are already present on
your local computer or samples that are stored in the `NCBI SRA
repository <https://www.ncbi.nlm.nih.gov/sra>`__.

-  For local sample files, identify a `glob
   pattern <https://en.wikipedia.org/wiki/Glob_(programming)>`__ that
   finds these files.
-  For samples on NCBI, identify the NCBI SRA run IDs of the samples you
   want to analyze. The Run IDs typically start with an SRR, ERR or DRR
   prefix. hese sample run IDs must be placed, one per line, in a file
   and the filename should be set in the ``remote_list_path`` of
   ``nextflow.config``.

.. note::

  The SRA Toolkit caches SRA files in your home directory by
  default. For large experiments this cache can become quite large, which
  may become an issue on some HPC clusters where each user is given a disk
  quota for their home directory. You can change the location of the cache
  by running ``vdb-config -i`` (see `SRA Toolkit
  Configuration <https://github.com/ncbi/sra-tools/wiki/Toolkit-Configuration>`__).

Step 2) Acquire Reference Genome Files
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Download the genome annotation/reference files. These will differ depending on
your organism that you are performing the analysis on. The files neccesary will
also be different  depending on if you decide to run Hisat2, Salmon or Kallisto.
The following instructions will go over how to make the proper files for each
of these programs.

For Hisat2
==========

You must aquire these 2 files from the internet:

-  A FASTA file containing the full genomic sequence (either
   pseudomolecules or scaffolds). Note, if your genome file is extremely
   large with hundreds of thousands of contigs/scaffolds, you may want
   to reduce the size of the FASTA file to contain only those
   contigs/scaffolds with predicted annotated genes. This will be converted into
   the hisat2 index files using the commands below.
-  A `GTF <https://uswest.ensembl.org/info/website/upload/gff.html>`__
   file containing the gene models.

   .. note::

     Sometimes a genome assembly does not
     provide a GTF file, but rather provides a
     `GFF3 <https://uswest.ensembl.org/info/website/upload/gff.html>`__
     file. You can convert the GFF file to a GTF file using the
     ``gffread`` program of
     `cufflinks <http://cole-trapnell-lab.github.io/cufflinks/file_formats/>`__,
     which you may have to download separately. An example command-line to
     convert a GFF3 to GTF is

     .. code:: bash

      gffread [gff_file] -T -o [gtf_file]

     where ``[gff_file]`` and ``[gtf_file]`` should be substituted for the
     names of your GFF3 and desired GTF file respectively.

Additional considerations:

-  You must have hisat2 index files of your genome sequence. These are
   constructed by using the ``hisat2-build`` command.

   .. code:: bash

    hisat2-build -f YourFastaFile.fna YourPrefix | tee > hisat2-build.log

-  The GTF file and the hisat2 index files must have the same prefix and
   this prefix must be identified in ``nextflow.config`` using the
   ``prefix`` parameter for ``hisat2-build``.

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


-  All of the genome annotation files must be in a directory and this
   directory must be identified in ``nextflow.config`` using the
   ``ref > path`` paramter.

For Salmon and Kallisto
=======================

The process for both of these is rather similar.

You must aquire this file from the internet:

- FASTA file containing the individual transcripts of the organism of interest.
  This should be formatted so that the name of each transcript is before the
  sequence of its corresponding transcipt. This is different from the FASTA file
  used for HISAT2.

This file then needs to be converted into the proper index file for that program
using the following commands. (assume the name of the transcript FASTA file is
``YourFasta.transcrits.fna``)

For Salmon:

.. code:: bash

  salmon index -t YourFasta.transcrits.fna -i YourFasta.transcrits.Salmon.indexed

For Kallisto:

.. code:: bash

  kallisto index -i YourFasta.transcripts.Kalisto.indexed YourFasta.transcrits.fna

.. note::

  If you are running GEMmaker with Docker Images, you will have to run these commands
  within the docker image. To get inside the docker image, NEEDS TO BE COMPLETED,
  CONTACT GEMmaker creators for details!


The GEMmaker repo contains an ``examples`` directory which contains
several small example setups. The ``LocalRunExample`` contains a
``reference`` directory whose files have the same prefix of ``CORG``.
This prefix is set for the ``prefix`` parameter in
``nextflow.config.example``. The ``RemoteRunExample`` also contains an
``SRA_IDS.txt`` file which contains a list of SRA fastq\_run\_IDs to
download from NCBI.

Step 3) Executing the Workflow
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

To execute the workflow on a local machine:

.. code:: bash

    nextflow run main.nf -profile standard

To resume a workflow in the event of a failure:

.. code:: bash

    nextflow run main.nf -profile standard -resume

To execute the workflow and generate trace, timeline and execution
reports:

.. code:: bash

    nextflow run main.nf -profile standard -with-report -with-timeline -with-trace

To execute the workflow on an HPC system you must edit
``nextflow.config`` and add an appropriate profile for your system.
Refer to the `Nextflow
documentation <https://www.nextflow.io/docs/latest/config.html#config-profiles>`__.
You can then use any of the above commands by changing the ``-profile``
argument to use your profile.

Performance Considerations
~~~~~~~~~~~~~~~~~~~~~~~~~~

For large experiments on an HPC system, it is important to make sure
that you are effectively utilizing the resources of the system. There
are a number of parameters in ``nextflow.config`` which can be used to
increase performance based on the capabilities of your system:

- ``params.execution.threads``: All processes which support multithreading
  (such as trimmomatic) will use this number of threads. This setting
  should be determined by the number of cores per node on your system; for
  example, if your system has nodes with 16 cores per node then you could
  set the number of threads to 16 to make full use of those nodes.

- ``params.execution.queue_size``: Nextflow will only run 100 processes at
  a time by default, but you may be able to increase this value based on
  the queue limits of your system.

Generating a Summary Report
~~~~~~~~~~~~~~~~~~~~~~~~~~~

The `MultiQC <http://multiqc.info>`__ tool will automatically generate a report
on how each process ran.

Generating the Gene Expression Matrix (GEM)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

After GEMmaker completes, the resulting GEM will be output to GEMmaker/output/GEM
(assuming that you didnt change the output directory location). This directory
contains the final GEM matrix, in raw, TPM and FPKM form. These can be used for
further analysis.

Using the GEM in other Analysis
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

DGE Analysis
============

The raw GEM matrix can be used for DGE analysis in edgeR and other DGE software.

Network Analysis
================

After construction of the GEM, network analysis can be performed.
`KINC <https://github.com/SystemsGenetics/KINC>`__ (Knowledge Independent
Network Construction) is a high performance gene co-expression  that can perform
Pearson's or Spearman's correlation with K-means or Gaussian mixture models. KINC
is a Qt/`ACE <https://github.com/SystemsGenetics/ACE>`__ application that is
capable of running on GPU's, making it fast and efficient.

.. |DOI| image:: https://zenodo.org/badge/114067776.svg
   :target: https://zenodo.org/badge/latestdoi/114067776
