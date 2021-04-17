Step 1: Prepare Genome Data
---------------------------

GEMmaker supports use of `Hisat2 <https://ccb.jhu.edu/software/hisat2/index.shtml>`_, `Kallisto <https://pachterlab.github.io/kallisto/>`_ and `Salmon <https://combine-lab.github.io/salmon/>`_, and allows you to select one of these tools to use for quantification of gene expression.  Each tool requires that transcript sequences of the genome are indexed prior to usage.

First, you must obtain the appropriate genome reference files and have them available on your local machine for indexing. Once you have obtained the files and placed them in a directory, you can index the genome by following the instructions in the sections below.

**You only need to index the files for the tool (i.e. kallisto, salmon or hisat2) that you would like GEMmaker to use.**


.. note::

  Sometimes the genome assembly for a species may have successive releases, with each improving either on the genome assembly or the genome structural annotations (i.e. identified genes) or both.  However, sometimes the functional annotation of genes may be lacking in most recent version as research communities developing the genome may release the genome for use prior to fully annotating it.  When you select a genome reference for use with GEMmaker, choose the assembly with the best functional annotations for genes, or ensure that functional annotations for the release you choose are available.  This will not affect the performance of GEMmaker but will affect downstream analyses such as identifying function of deferentially expressed genes or modules of genes from a co-expression network.

The instructions below provide examples for indexing the *Arabidopsis thaliana* genome as if it were obtained from `Ensemble Plants <http://plants.ensembl.org/>`_.

Kallisto
''''''''
Kallisto is the default tool that GEMmaker uses for gene transcript quantification. For Kallisto, you need the nucleotide sequences of all transcripts for your species in FASTA format.  Kallisto recommends, for example, using the cDNA FASTA files similar to what you find on `Ensembl genomes <http://ensemblgenomes.org/>`__. After obtaining the file, you must use the ``kallisto index`` command to index the file. For more usage information please, see the `Kallisto Online Manual <https://pachterlab.github.io/kallisto/manual>`_.

For example, to retrieve the Arabidopsis cDNA file:

.. code-block:: bash

  wget ftp://ftp.ensemblgenomes.org/pub/plants/release-50/fasta/arabidopsis_thaliana/cdna/Arabidopsis_thaliana.TAIR10.cdna.all.fa.gz


Index Kallisto using Singularity
................................
If you do not have Kallisto indexes already prepared for your reference genome, you can use the GEMmaker docker image to perform the indexing. For example, you can use Singularity in the following way:

.. code-block:: bash

  singularity exec -B ${PWD} docker://systemsgenetics/gemmaker kallisto index -i Arabidopsis_thaliana.TAIR10.kallisto.indexed Arabidopsis_thaliana.TAIR10.cdna.all.fa.gz

The command above uses the ``gemmaker/gemmaker`` image that was built for GEMmaker.  The image will be downloaded if it does not already exist on your machine.  The command above uses the ``-B ${PWD}`` argument to automatically mount the current directory onto the same directory in the image. From there the Kallisto index command can be executed.

Index Kallisto using Docker
...........................
If you do not have Kallisto indexes already prepared for your reference genome, you can use the GEMmaker docker image to perform the indexing. For example, you can use Docker in the following way:

.. code-block:: bash

  docker run -v ${PWD}:/reference -u $(id -u ${USER}):$(id -g ${USER}) systemsgenetics/gemmaker /bin/bash -c "cd reference; kallisto index -i Arabidopsis_thaliana.TAIR10.kallisto.indexed Arabidopsis_thaliana.TAIR10.cdna.all.fa.gz"

The command above uses the ``gemmaker/gemmaker`` image that was built for GEMmaker.  The image will be downloaded if it does not already exist on your machine.  The ``-v ${PWD}:/references`` argument instructs Docker to mount the current directory (i.e.: ``${PWD}``) onto a new directory in the image named ``/reference`` and gives the image access to the transcript file for indexing. From there the Kallisto index command can be executed.  The ``-c`` argument provides the Kallisto index command needed to index the files.  The ``-u $(id -u ${USER}):$(id -g ${USER})`` argument instructs Docker to run Kallisto as you rather than the system root user.

Salmon
''''''
Salmon is similar in terms of performance and results as Kallisto. For Salmon, you need the nucleotide sequences of all transcripts for your species in FASTA format.  Be sure to find a FASTA file containing cDNA sequences. After obtaining the file, you must use the ``salmon index`` command to index the file. For more usage information please, see the `Salmon Online Manual <https://salmon.readthedocs.io/en/latest/index.html>`_.

As an example, to retrieve the Arabidopsis cDNA file:

.. code-block:: bash

  wget ftp://ftp.ensemblgenomes.org/pub/plants/release-50/fasta/arabidopsis_thaliana/cdna/Arabidopsis_thaliana.TAIR10.cdna.all.fa.gz


Index Salmon using Singularity
..............................
If you do not have Salmon indexes already prepared for your reference genome, you can use the GEMmaker docker image to perform the indexing. For example, you can use Singularity in the following way:

.. code-block:: bash

   singularity exec -B ${PWD} docker://systemsgenetics/gemmaker salmon index index -t Arabidopsis_thaliana.TAIR10.cdna.all.fa.gz -i Arabidopsis_thaliana.TAIR10.salmon.indexed


The command above uses the ``systemsgenetics/gemmaker`` image to index the transcripts.  The image will be downloaded if it does not already exist on your machine.  The command above uses the ``-B ${PWD}`` argument to automatically mount the current directory onto the same directory in the image. From there the Salmon index command can be executed.

Index Salmon using Docker
.........................
If you do not have Salmon indexes already prepared for your reference genome, you can use the GEMmaker docker image to perform the indexing. For example, you can use Docker in the following way:

.. code-block:: bash

  docker run -v ${PWD}:/reference -u $(id -u ${USER}):$(id -g ${USER}) systemsgenetics/gemmaker /bin/bash -c "cd /reference; salmon index index -t Arabidopsis_thaliana.TAIR10.cdna.all.fa.gz -i Arabidopsis_thaliana.TAIR10.salmon.indexed"

The command above uses the ``systemsgenetics/gemmaker`` image that was built by the GEMmaker development team to index the transcripts.  The image will be downloaded if it does not already exist on your machine.   The ``-v ${PWD}:/references`` argument instructs Docker to mount the current directory (i.e.: ``${PWD}``) onto a new directory in the image named ``/reference`` and gives the image access to the transcript file for indexing.  The ``-c`` argument provides the Salmon index command needed to index the files.  The ``-u $(id -u ${USER}):$(id -g ${USER})`` argument instructs Docker to run Salmon as you rather than the system root user.

Hisat2
''''''
Hisat2 is different from Kallisto and Salmon in that it requires multiple steps that include alignment of RNA-seq reads to a genomic reference sequence followed by quantification of expression using the tool `StringTie <https://ccb.jhu.edu/software/stringtie/>`__. You must therefore obtain the following files:

-  A FASTA file containing the full genomic sequence in FASTA format (either pseudomolecules or scaffolds).
-  A `GTF <https://uswest.ensembl.org/info/website/upload/gff.html>`__ file containing the gene models.

As an example, to retreive the Arabidopsis files:

.. code-block:: bash

  wget ftp://ftp.ensemblgenomes.org/pub/plants/release-50/fasta/arabidopsis_thaliana/dna/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.gz
  gunzip Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.gz

  wget ftp://ftp.ensemblgenomes.org/pub/plants/release-50/gff3/arabidopsis_thaliana/Arabidopsis_thaliana.TAIR10.50.gff3.gz
  gunzip Arabidopsis_thaliana.TAIR10.50.gff3.gz

.. note::
  If your genome file is extremely large with hundreds of thousands of contigs/scaffolds, you may want to reduce the size of the FASTA file to contain only those contigs/scaffolds with predicted annotated genes.

Sometimes a genome assembly does not provide a GTF file, but rather provides a `GFF3 <https://uswest.ensembl.org/info/website/upload/gff.html>`__ file. This is the case for the Arabidopsis genome provided by Ensemble You can convert the GFF file to a GTF file using the `gffread <https://github.com/gpertea/gffread>`__.  Examples for using gffread are provdied below.


Index Hisat2 using Singularity
..............................
If you do not have a GTF or Hisat2 indexes already prepared for your reference genome, you can use the GEMmaker docker image to create the GTF and perform the indexing. For example, you can use Singularity in the following way:

To create the GTF file:

.. code-block:: bash

  singularity exec -B ${PWD} docker://systemsgenetics/gemmaker  gffread Arabidopsis_thaliana.TAIR10.50.gff3.gz -T -o Arabidopsis_thaliana.TAIR10.gtf

To index the reference:

.. code-block:: bash

   singularity exec -B ${PWD} docker://systemsgenetics/gemmaker hisat2-build -f Arabidopsis_thaliana.TAIR10.dna.toplevel.fa Arabidopsis_thaliana.TAIR10

The following describes the meaning of the arguments in the command-line above:

The command above uses the ``systemsgenetics/gemmaker`` image.  The image will be downloaded if it does not already exist on your machine.  The command above uses the ``-B ${PWD}`` argument to automatically mount the current directory onto the same directory in the image. From there the Hisat2 index command can be executed.

Index Hisat2 using Docker
.........................
If you do not have a GTF or Hisat2 indexes already prepared for your reference genome, you can use the GEMmaker docker image to create the GTF and perform the indexing. For example, you can use Docker in the following way:


To create the GTF file:

.. code-block:: bash

  docker run -v ${PWD}:/reference -u $(id -u ${USER}):$(id -g ${USER}) systemsgenetics/gemmaker /bin/bash -c "cd /reference; gffread Arabidopsis_thaliana.TAIR10.50.gff3 -T -o Arabidopsis_thaliana.TAIR10.gtf"

To index the reference:

.. code-block:: bash

  docker run -v ${PWD}:/reference -u $(id -u ${USER}):$(id -g ${USER}) systemsgenetics/gemmaker  /bin/bash -c "cd /reference; hisat2-build -f Arabidopsis_thaliana.TAIR10.dna.toplevel.fa Arabidopsis_thaliana.TAIR10"

The command above uses the ``systemsgenetics/gemmaker`` image.  The image will be downloaded if it does not already exist on your machine.   The ``-v ${PWD}:/reference`` argument instructs Docker to mount the current directory (i.e.: ``${PWD}``) onto a new directory in the image named ``/references`` and gives the image access to the transcript file for indexing.  The ``-c`` argument provides the Salmon index command needed to index the files.  The ``-u $(id -u ${USER}):$(id -g ${USER})`` argument instructs Docker to run ``hisat2-build`` as you rather than the system root user.
