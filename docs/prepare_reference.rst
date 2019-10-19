Step 2: Prepare Genome Data
---------------------------

.. warning::

  Before proceeding please remember to clear out the ``input`` directory of example data as described on the :doc:`./usage` page, or you may accidentally include example data with your own data.

GEMmaker supports use of `Hisat2 <https://ccb.jhu.edu/software/hisat2/index.shtml>`_, `Kallisto <https://pachterlab.github.io/kallisto/>`_ and `Salmon <https://combine-lab.github.io/salmon/>`_, and allows you to select one of these tools to use for quantification of gene expression.  Each tool requires that transcript sequences of the genome are indexed prior to usage.

First, you must obtain the appropriate genome reference files and have them available on your local machine for indexing. By default, GEMmaker expects to find all genome reference files in the ``input/references`` directory of GEMmaker. Once you have obtained the files and placed them in the ``input/references`` directory, you can index the genome by following the instructions in the sections below. **You only need to index the files for the tool you would like GEMmaker to use.**


.. note::

  Sometimes the genome assembly for a species may have successive releases, with each improving either on the genome assembly or the genome structural annotations (i.e. identified genes) or both.  However, sometimes the functional annotation of genes may be lacking in most recent version as research communities developing the genome may release the genome for use prior to fully annotating it.  When you select a genome reference for use with GEMmaker, choose the assembly with the best functional annotations for genes, or ensure that functional annotations for the release you choose are available.  This will not affect the performance of GEMmaker but will affect downstream analyses such identifying function of deferentially expressed genes or modules of genes from a co-expression network.


Kallisto
''''''''
Kallisto is the default tool that GEMmaker uses for gene transcript quantification. For Kallisto, you need the nucleotide sequences of all transcripts for your species in FASTA format.  Kallisto recommends, for example, using the cDNA FASTA files similar to what you find on `Ensembl genomes <http://ensemblgenomes.org/>`__. After obtaining the file, you must use the ``kallisto index`` command to index the file. For more usage information please, see the `Kallisto Online Manual <https://pachterlab.github.io/kallisto/manual>`_.

The instructions below provide examples for indexing the fake CORG example data that comes with GEMmaker. These data are found in the ``input/references`` directory of GEMmaker.  Here will be shown how to index the ``CORG.transcripts`` FASTA file into an index file named ``CORG.transcripts.Kallisto.indexed``.


If Kallisto is installed locally
................................

If you have Kallisto installed locally, you can create the index, using the following command:

.. code-block:: bash

   kallisto index -i CORG.transcripts.Kallisto.indexed CORG.transcripts


Index Kallisto using Singularity
................................
If you do not have Kallisto installed locally, but you have Singularity installed, you can use the GEMmaker Kallisto docker image to perform the indexing.

.. code-block:: bash

  singularity exec -B ${PWD} docker://gemmaker/kallisto:0.45.0-1.1 kallisto index -i CORG.transcripts.Kallisto.indexed CORG.transcripts

The following describes the meaning of the arguments in the command-line above:

- `singularity exec -B ${PWD} docker://gemmaker/kallisto:0.45.0-1.1 kallisto  index`: Downloads the kallisto singularity image to the current directory and runs it
- `-i CORG.transcripts.Kallisto.indexed`: Name of output (indexed) Kallisto File
- `CORG.transcripts`: Name of input transcript files

The command above uses the ``gemmaker/kallisto:0.45.0-1.1`` image that was built by the GEMmaker development team to index the transcripts.  The image will be downloaded if it does not already exist on your machine.  The command above uses the ``-B ${PWD}`` argument to automatically mount the current directory onto the same directory in the image. From there the Kallisto index command can be executed.

Index Kallisto using Docker
...........................
If you do not have Kallisto installed locally, but you have Docker installed, you can use the GEMmaker Kallisto docker image to perform the indexing.

.. code-block:: bash

  docker run -v ${PWD}:/GEMmaker/input/references -u $(id -u ${USER}):$(id -g ${USER}) gemmaker/kallisto:0.45.0-1.1 /bin/bash -c "cd /GEMmaker/input/references; kallisto index -i CORG.transcripts.Kallisto.indexed CORG.transcripts"

The command above uses the ``gemmaker/kallisto:0.45.0-1.1`` image that was built by the GEMmaker development team to index the transcripts.  The image will be downloaded if it does not already exist on your machine.   The ``-v ${PWD}:/GEMmaker/input/references`` argument instructs Docker to mount the current directory (i.e.: ``${PWD}``) onto a new directory in the image named ``/GEMMaker/input/references`` and gives the image access to the transcript file for indexing.  The ``-c`` argument provides the Kallisto index command needed to index the files.  The ``-u $(id -u ${USER}):$(id -g ${USER})`` argument instructs Docker to run Kallisto as you rather than the system root user.

Salmon
''''''
Salmon is similar in terms of performance and results as Kallisto. For Salmon, you need the nucleotide sequences of all transcripts for your species in FASTA format.  Be sure to find a FASTA file containing cDNA sequences. After obtaining the file, you must use the ``salmon index`` command to index the file. For more usage information please, see the `Salmon Online Manual <https://salmon.readthedocs.io/en/latest/index.html>`_.

The instructions below provide examples for indexing the fake CORG example data that comes with GEMmaker. These data are found in the ``input/references`` directory of GEMmaker.  Here will be shown how to index the ``CORG.transcripts`` FASTA file into an index file named ``CORG.transcripts.Kallisto.indexed``.


If Salmon is installed locally
..............................

If you have Kallisto installed locally, you can create the index, using the following command:

.. code-block:: bash

  salmon index -t CORG.transcripts -i CORG.transcripts.Salmon.indexed

Index Salmon using Singularity
..............................
If you do not have Salmon installed locally, but you have Singularity installed, you can use the GEMmaker Salmon docker image to perform the indexing.

.. code-block:: bash

   singularity exec -B ${PWD} docker://gemmaker/salmon:0.12.0-1.1 salmon index -t CORG.transcripts -i CORG.transcripts.Salmon.indexed

The following describes the meaning of the arguments in the command-line above:

- `singularity exec -B ${PWD} docker://gemmaker/salmon:0.12.0-1.1 salmon index`: Downloads the salmon singularity image to the current directory and runs it
- `-t CORG.transcripts` : Name of input transcript files
- `-i CORG.transcripts.Salmon.indexed` Name of output (indexed) Salmon fiel


The command above uses the ``gemmaker/salmon:0.12.0-1.1`` image that was built by the GEMmaker development team to index the transcripts.  The image will be downloaded if it does not already exist on your machine.  The command above uses the ``-B ${PWD}`` argument to automatically mount the current directory onto the same directory in the image. From there the Salmon index command can be executed.

Index Salmon using Docker
.........................
If you do not have Salmon installed locally, but you have Docker installed, you can use the GEMmaker Salmon docker image to perform the indexing.

.. code-block:: bash

  docker run -v ${PWD}:/GEMmaker/input/references -u $(id -u ${USER}):$(id -g ${USER}) gemmaker/salmon:0.12.0-1.1 /bin/bash -c "cd /GEMmaker/input/references; salmon index -t CORG.transcripts -i CORG.transcripts.Salmon.indexed"

The command above uses the ``gemmaker/salmon:0.12.0-1.1`` image that was built by the GEMmaker development team to index the transcripts.  The image will be downloaded if it does not already exist on your machine.   The ``-v ${PWD}:/GEMmaker/input/references`` argument instructs Docker to mount the current directory (i.e.: ``${PWD}``) onto a new directory in the image named ``/GEMMaker/input/references`` and gives the image access to the transcript file for indexing.  The ``-c`` argument provides the Salmon index command needed to index the files.  The ``-u $(id -u ${USER}):$(id -g ${USER})`` argument instructs Docker to run Salmon as you rather than the system root user.

Hisat2
''''''
Hisat2 is different from Kallisto and Salmon in that it requires multiple steps that include alignment of RNA-seq reads to a genomic reference sequence followed by quantification of expression using the tool `StringTie <https://ccb.jhu.edu/software/stringtie/>`__. You must therefore obtain the following files:

-  A FASTA file containing the full genomic sequence in FASTA format (either pseudomolecules or scaffolds).
-  A `GTF <https://uswest.ensembl.org/info/website/upload/gff.html>`__ file containing the gene models.

.. note::
  If your genome file is extremely large with hundreds of thousands of contigs/scaffolds, you may want to reduce the size of the FASTA file to contain only those contigs/scaffolds with predicted annotated genes.

.. note::
  Sometimes a genome assembly does not provide a GTF file, but rather provides a `GFF3 <https://uswest.ensembl.org/info/website/upload/gff.html>`__ file. You can convert the GFF file to a GTF file using the ``gffread`` tool from `cufflinks <http://cole-trapnell-lab.github.io/cufflinks/file_formats/>`__, which you may have to download separately. Here is an example command-line to convert a GFF3 to GTF:

  .. code:: bash

    gffread <gff_file> -T -o <gtf_file>

  The arguments ``<gff_file>`` and ``<gtf_file>`` should be substituted for the names of your GFF3 and desired GTF file respectively.

After obtaining your genome sequence file it must be indexed. These are constructed by using the ``hisat2-build`` command.   Hisat indexes are contained in multiple files with the same prefix and a ``.ht2`` extension.  The following provides instructions for indexing the genome using hisat2.


If Hisat2 is installed locally
..............................
If Hisat2 is installed locally, you can create the indexes, using the following command.

  .. code:: bash

    hisat2-build -f CORG.fna CORG


Index Hisat2 using Singularity
..............................
If you do not have Hisat2 installed locally, but you have Singularity installed, you can use the GEMmaker Hisat2 docker image to perform the indexing.

.. code-block:: bash

   singularity exec -B ${PWD} docker://gemmaker/hisat2:2.1.0-1.1 hisat2-build -f CORG.fna CORG

The following describes the meaning of the arguments in the command-line above:

- `singularity exec -B ${PWD} docker://gemmaker/hisat2:2.1.0-1.1 hisat2-build`: Downloads the hisat2 singularity image to the current directory and runs it
- `-f CORG.fna` : The input genome fasta file
- `CORG`: The name of the reference organism. This will be used as a prefix for all the reference files.

The command above uses the ``gemmaker/hisat2:2.1.0-1.1`` image that was built by the GEMmaker development team to index the transcripts.  The image will be downloaded if it does not already exist on your machine.  The command above uses the ``-B ${PWD}`` argument to automatically mount the current directory onto the same directory in the image. From there the Hisat2 index command can be executed.

Index Hisat2 using Docker
.........................
If you do not have Salmon installed locally, but you have Docker installed, you can use the GEMmaker Salmon docker image to perform the indexing.

.. code-block:: bash

  docker run -v ${PWD}:/GEMmaker/input/references -u $(id -u ${USER}):$(id -g ${USER}) gemmaker/hisat2:2.1.0-1.1 /bin/bash -c "cd /GEMmaker/input/references; hisat2-build -f CORG.fna CORG"

The command above uses the ``gemmaker/hisat2:2.1.0-1.1`` image that was built by the GEMmaker development team to index the transcripts.  The image will be downloaded if it does not already exist on your machine.   The ``-v ${PWD}:/GEMmaker/input/references`` argument instructs Docker to mount the current directory (i.e.: ``${PWD}``) onto a new directory in the image named ``/GEMMaker/input/references`` and gives the image access to the transcript file for indexing.  The ``-c`` argument provides the Salmon index command needed to index the files.  The ``-u $(id -u ${USER}):$(id -g ${USER})`` argument instructs Docker to run ``hisat2-build`` as you rather than the system root user.
