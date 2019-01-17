.. _running_the_examples:

Testing GEMmaker
----------------

The GEMmaker example run will with a small set of local files (provided in the GEMmaker repo) and a remote file (automatically downloaded from `NCBI's SRA repository <https://www.ncbi.nlm.nih.gov/sra>`__). These samples are extremely small and are only meant to demonstrate usage for a mixed set of local and remote files.

Once the :ref:`software_prerequisites` are installed, you can run GEMmaker with the example data. If you have not already done so, you must create your own config file from the example:

.. code:: bash
  cp nextflow.config.example nextflow.config

The example ``nextflow.config`` is already configured with the proper settings to run the example described above. To execute the workflow run the following:

.. code:: bash
  nextflow run main.nf -profile standard,localDocker

You should see an output that looks like this: :ref:`example_output`

Explanation of the Example Run
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This example uses the imaginary organism "Cool Organism" (CORG). For the local example, we use a set of 3 artificially made RNA-seq runs made for CORG. CORG has a very small "genome" of only 2,336 nucleotides, 3 "chromosomes" and 6 "genes". The 6 genes are named ``gene_Alpha``, ``gene_Beta``, ``gene_Zeta``, ``gene_Gamma``, ``gene_Delta``, ``gene_Epsilon``.

For the remote example, GEMmaker automatically downloads a very small RNA-seq file from NCBI. This dataset is from an uncharacterized bacteria, but luckly, CORG shares 3 of the genes with this bacteria so we can use CORG's reference file (pretend that the remote file is also for CORG, we are just using it becasue it is an unusually small file, which makes it an ideal example).

Using Salmon or Kallisto
~~~~~~~~~~~~~~~~~~~~~~~~

The example uses Hisat2 by default. If you would like to use Salmon or Kallisto instead you must edit ``nextflow.config`` and change the alignment type. In the GEMmaker directory, edit ``nextflow.config`` using your favorite text editor. Here we use `vim <https://www.vim.org/>`__ on the command line:

.. code:: bash
  vim nextflow.config

Then edit ``params.software.alignment`` in the config file. Change to ``1`` for Kallisto, and ``2`` for salmon. For example, to use Kallisto:

.. code:: bash
  //
  // hisat2 = 0
  // Kallisto = 1
  // Salmon = 2
  //
  alignment = 1

After that, save your file and run the worklow:

.. code:: bash
  nextflow run main.nf -profile standard,localDocker

Explanation of the Inputs
~~~~~~~~~~~~~~~~~~~~~~~~~

The inputs for the example run are in the ``examples`` directory, and they consist of the reference directory and two data directories for local and remote samples.

Reference directory
===================

The reference directory for the example is located at:

.. code:: bash
  GEMmaker/examples/reference/

This directory contains the

- made up reference genome file (``CORG.fna``),
- `GTF <https://uswest.ensembl.org/info/website/upload/gff.html>`__ file (``CORG.gtf``)
- hisat index files (``CORG.?/ht2``).
- kallisto index file (``CORG.transcripts.Kallisto.indexed``)
- salmon index directory (``CORG.transcripts.Salmon.indexed/``)
- ``COMMANDS.sh`` explaining how each of these were generated

These are the files needed to run hisat2, kallisto, and salmon on the CORG data.

Data directories
================

There are two sample data directories:

For local runs:

.. code:: bash
  GEMmaker/examples/LocalRunExample/

For remote runs:

.. code:: bash
  GEMmaker/examples/RemoteRunExample/

The `LocalRunExample` directory contains three `FASTQ <https://en.wikipedia.org/wiki/FASTQ_format>`__ files for CORG containing RNA-seq data. These are examples of local unpaired data, and are each in a directory of their own. The file naming format for these reads is "?\_1.fastq" where the "?" is the number of the sample. GEM-maker finds these files through the glob pattern assigned to the ``local_samples_path`` in ``nextflow.config``.

The `RemoteRunExample` directory contains the file ``SRA_IDs.txt`` which contains a list of names for remote files to be downloaded by GEMmaker from `NCBI's SRA repository <https://www.ncbi.nlm.nih.gov/sra>`__. In this case, there is only one run ID.

Explanation of the Outputs
~~~~~~~~~~~~~~~~~~~~~~~~~~

Once executed, the example should create a directory called ``output`` with several sub-directories. Four of these directories correspond to each sample (3 local, 1 remote), and each of these contains the files generated for that sample. In each sample directory you will find the following files:

- ``fastq``: The fastq reads file for the experiment.
- ``fastqc``: 6 or 12 files (depending on paired or unpaired data) from fastqc. FastQC is configured to check files before and after trimmomatic.
- ``bam``: Binary alignment file.
- ``ga``: Expression level transcript abundance.
- ``fpkm``: Two-column version of the ``ga`` file with only gene and FPKM value.
- ``tpm``: Two-column version of the ``ga`` file with only gene and TPM value.

The other directories are the ``reports`` directory and the ``GEM`` directory. The ``reports`` directory will contain a ``multiqc_report.html`` file that provides several statistics about the run.

.. figure:: /images/MultiQC_Report.png
  :alt: MultiQC_Report

Figure 1: Image of the start of the report for the example run when run with Hisat2.

The ``GEMs`` directory contains the final gene-expression matrices (GEMs), in raw, TPM and FPKM form. These can be used as input to other analyses.

The output of GEMmaker can be used for several different analyses. The FPKM files can be combined into an expression matrix and then visualized using a heatmap. The heatmap below consists of the FPKM values from the local examples divided by 1000. We can see that ``gene_Zeta`` remained constant across all three samples, ``gene_Beta`` decreased, and ``gene_Alpha`` increased.

.. figure:: /images/heatmap.png
  :alt: heatmap

Figure 2: Heatmap of FPKM values from local samples.
