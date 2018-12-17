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

Step 2) Aquire Reference Genome Files
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Download the genome annotation/reference files. These will differ depending on
your organism that you are performing the analysis on. You must have the
following:

-  A FASTA file containing the full genomic sequence (either
   pseudomolecules or scaffolds). Note, if your genome file is extremely
   large with hundreds of thousands of contigs/scaffolds, you may want
   to reduce the size of the FASTA file to contain only those
   contigs/scaffolds with predicted annotated genes.
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

The `MultiQC <http://multiqc.info>`__ tool can be used with GEMmaker to
generate a summary report of results from Trimmomatic, Hisat2 and
samtools. This report allows you to explore the quality of the data,
trimming and alignments. To generate the report you must have `MultiQC
installed <http://multiqc.info/docs/#installing-multiqc>`__. Once
installed, you can generate the report with the following command inside
of the GEMmaker directory where your workflow was executed:

.. code:: bash

    multiqc .

Generating the Gene Expression Matrix (GEM)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

After GEMmaker completes, the results for each sample are stored in a
directory specific to that sample. The final output for each sample is a
Gene Expression Vector (GEV) in the form of an FPKM or TPM file. To
compile all GEVs into a Gene Expression Matrix (GEM) you can use the
``create_GEM.py`` script in the ``scripts`` directory.

To see help documentation for this script:

.. code:: bash

    python ./scripts/create_GEM.py -h

To create a GEM file from the TPM files produced by GEMmaker:

.. code:: bash

    python ./scripts/create_GEM.py --source ./ --type TPM --prefix my_project

The script will produce a GEM file called ``my_project.GEM.TPM.txt``. Be
sure to change ``my_project`` to a meaningful prefix for your project.

You can combine the results of multiple GEMmaker runs into a single GEM
by providing a list of directories to the ``--source`` argument. This
feature may be useful if you split a set of input files into several
GEMmaker runs and now you need to combine then. The script will produce
a file named ``GEM.txt`` in the working directory.

Using the GEM in other Analysis
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

DGE Analysis
============
Need to do research on:

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
