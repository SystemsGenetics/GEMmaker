Data Preparation
----------------
Prior to using GEMmaker you must prepare intput data. This includes locating and indexing the reference genome files and identifying the input RNA-seq data files you wish to process.

Preparing Genome Reference Files
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
GEMmaker supports use of `Hisat2 <https://ccb.jhu.edu/software/hisat2/index.shtml>`_, `Kallisto <https://pachterlab.github.io/kallisto/>`_ and `Salmon <https://combine-lab.github.io/salmon/>`_, and allows you to select one of these tools to use for quantification of gene expression.  Each tool requires that transcript sequences of the genome are indexed prior to usage.

First, you must obtain the appropriate genome reference files and have them available on your local machine for indexing.  Then, you can index the genome by following the instructions in the sections below.  You only need to index the files for the tool you would like GEMmaker to use. For these instructions you must have the tool downloaded and installed on your local machine, or you must have Docker.


.. note::
  Sometimes the genome assembly for a species may have successive releases, with each improving either on the genome assembly or the genome structural annotations (i.e. identified genes) or both.  However, sometimes the functional annotation of genes may be lacking in most recent version as research communities developing the genome may release the genome for use prior to fully annotating it.  When you select a genome reference for use with GEMmaker, choose the assembly with the best functional annotations for genes, or ensure that functional annotations for the release you choose are available.  This will not affect the performance of GEMmaker but will affect downstream analyses such identifying function of deferentially expressed genes or modules of genes from a co-expression network.

Hisat2
""""""

Salmon
""""""

Kallisto
""""""""
For Kallisto, you need the nucleotide sequences of all transcripts for your species in FASTA format.  Often, such a file accompanies the genome assembly. You must use the ``kallisto index`` command to index the file. For example, if the FASTA file is named `transcripts.fna`, then the following command would index this file for use with kallisto:

.. code-block:: bash

   kallisto index -i transcripts.indexed transcripts.fna

The file `transcripts.indexed` contains the indexes for teh sequences and can now be used with GEMmaker. For more usage information please, see the `Kallisto Online Manual <https://pachterlab.github.io/kallisto/manual>`_.

If you do not have Kallisto installed on your local machine, but you have Docker installed, you can use the GEMmaker Kallisto docker image to perform the indexing.  For a FASTA filed named `transcripts.fna`, the following command will perform the indexing as long as it is run within the same directory as the FASTA file:

.. code-block:: bash

  docker run -v ${PWD}:/mnt systemsgenetics/kallisto:0.45.0  kallisto index -i /mnt/transcripts.indexed /mnt/transcripts.fna

In the command above the ``-v`` option instructs Docker to mount the current directory (specirifed by ``${PWD}``) onto the ``/mnt`` directory in the image. This allows the image to see the FASTA file.
