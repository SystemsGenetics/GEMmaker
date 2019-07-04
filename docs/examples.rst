.. _examples:

Example
-------

The GEMmaker example consists of a small set of local files (provided in the GEMmaker repo) and a remote file (automatically downloaded from `NCBI's SRA repository <https://www.ncbi.nlm.nih.gov/sra>`__). These samples are extremely small and are only meant to demonstrate usage for a mixed set of local and remote files.

Once GEMmaker and its dependencies have been installed, you can run GEMmaker with the example data. The ``nextflow.config`` file is already configured to run this example. Execute the workflow using the following command:

.. code:: bash

  nextflow run main.nf -profile standard -with-singularity

You should see an output that looks like this:

.. code:: bash

  N E X T F L O W  ~  version 19.04.1
  Launching `main.nf` [spontaneous_aryabhata] - revision: 333c8c3e4f

  ===================================
   G E M M A K E R   P I P E L I N E
  ===================================

  Workflow Information:
  ---------------------
    Project Directory:  /home/usr/Documents/GEMmaker
    Launch Directory:   /home/usr/Documents/GEMmaker
    Work Directory:     /home/usr/Documents/GEMmaker/work
    Config Files:       [/home/usr/Documents/GEMmaker/nextflow.config]
    Container Engine:   singularity
    Profile(s):         standard


  Input Parameters:
  -----------------
    Remote fastq list path:     /home/usr/Documents/GEMmaker/examples/RemoteRunExample/SRA_IDs.txt
    Local sample glob:          /home/usr/Documents/GEMmaker/examples/LocalRunExample/Sample*/*_{1,2}.fastq
    Reference genome path:      /home/usr/Documents/GEMmaker/examples/reference/
    Reference genome prefix:    CORG


  Output Parameters:
  ------------------
  Output directory:           output
  Publish SRA:                false
  Publish downloaded FASTQ:   false
  Publish trimmed FASTQ:      false
  Publish BAM:                false
  Publish Gene Abundance:     false
  Publish GTF_GA:             false
  Publish RAW:                true
  Publish FPKM:               true
  Publish TPM:                true
  MultiQC:                    true
  Create GEM:                 true


  Execution Parameters:
  ---------------------
    Queue size:                 100


  Software Parameters:
  --------------------
    Trimmomatic clip path:      /home/usr/Documents/GEMmaker/files/fasta_adapter.txt
    Trimmomatic minimum ratio:  0.7

  [warm up] executor > local
  executor >  local (70)
  [c5/16ba32] process > write_stage_files      [100%] 4 of 4 ✔
  [b8/f56312] process > retrieve_sra_metadata  [100%] 1 of 1 ✔
  [6b/3645f9] process > start_first_batch      [100%] 1 of 1 ✔
  [b6/aca093] process > read_sample_file       [100%] 4 of 4 ✔
  [f6/cd627d] process > fastqc_1               [100%] 4 of 4 ✔
  [cd/bfe65b] process > prefetch               [100%] 1 of 1 ✔
  [0a/808a19] process > trimmomatic            [100%] 4 of 4 ✔
  [1b/e4d65d] process > fastq_dump             [100%] 1 of 1 ✔
  [a9/80eb92] process > fastqc_2               [100%] 4 of 4 ✔
  [62/45f1d2] process > hisat2                 [100%] 4 of 4 ✔
  [4e/88ae16] process > fastq_merge            [100%] 1 of 1 ✔
  [54/10abb0] process > clean_sra              [100%] 1 of 1 ✔
  [87/3ddce8] process > samtools_sort          [100%] 4 of 4 ✔
  [20/356849] process > clean_downloaded_fastq [100%] 1 of 1 ✔
  [c6/0bc5dc] process > samtools_index         [100%] 4 of 4 ✔
  [83/c51c00] process > clean_sam              [100%] 4 of 4 ✔
  [69/2f1adb] process > clean_trimmed_fastq    [100%] 4 of 4 ✔
  [0e/077fe4] process > stringtie              [100%] 4 of 4 ✔
  [3e/ed9160] process > hisat2_fpkm_tpm        [100%] 4 of 4 ✔
  [2a/5a840b] process > clean_bam              [100%] 4 of 4 ✔
  [ca/91419b] process > next_sample            [100%] 4 of 4 ✔
  [ba/db160c] process > clean_stringtie_ga     [100%] 4 of 4 ✔
  [b2/3ac691] process > multiqc                [100%] 1 of 1 ✔
  [a0/6aeb80] process > create_gem             [100%] 1 of 1 ✔
  [81/a3f835] process > clean_merged_fastq     [100%] 1 of 1 ✔
  Completed at: 02-Jul-2019 13:20:08
  Duration    : 47.9s
  CPU hours   : (a few seconds)
  Succeeded   : 70




Additionally, you should see a directory called ``output`` with the following subdirectories:

.. code:: bash

  output/
    1/
    2/
    3/
    GEMs/
    reports/
    SRX218012/

The "CORG" Example
~~~~~~~~~~~~~~~~~~

This example uses the imaginary organism "Cool Organism" (CORG). For the local example, we use a set of 3 artificially made RNA-seq runs made for CORG. CORG has a very small "genome" of only 2,336 nucleotides, 3 "chromosomes" and 6 "genes". The 6 genes are named ``gene_Alpha``, ``gene_Beta``, ``gene_Zeta``, ``gene_Gamma``, ``gene_Delta``, ``gene_Epsilon``.

For the remote example, GEMmaker automatically downloads a very small RNA-seq file from NCBI. This dataset is from an uncharacterized bacteria, but luckily, CORG shares 3 of the genes with this bacteria so we can use CORG's reference file (pretend that the remote file is also for CORG, we are just using it becasue it is an unusually small file, which makes it an ideal example).

Using Salmon or Kallisto
~~~~~~~~~~~~~~~~~~~~~~~~

The example uses Hisat2 by default. If you would like to use Salmon or Kallisto instead, you must edit ``nextflow.config`` and change the alignment type. In the GEMmaker directory, edit ``nextflow.config`` using your favorite text editor. Here we use `vim <https://www.vim.org/>`__ on the command line:

.. code:: bash

  vim nextflow.config

Then edit ``params.software.alignment`` in the config file. Change to ``1`` for Kallisto, and ``2`` for Salmon. For example, to use Kallisto:

.. code:: bash

  //
  // hisat2 = 0
  // Kallisto = 1
  // Salmon = 2
  //
  alignment = 1

Then save your file and run the worklow:

.. code:: bash

  nextflow run main.nf -profile standard -with-docker

Explanation of the Inputs
~~~~~~~~~~~~~~~~~~~~~~~~~

The inputs for the example run are in the ``examples`` directory, and they consist of the reference directory and two data directories for local and remote samples.

Reference directory
===================

The reference directory for the example is located at:

.. code:: bash

  GEMmaker/examples/reference/

This directory contains the

- reference genome file (``CORG.fna``),
- `GTF <https://uswest.ensembl.org/info/website/upload/gff.html>`__ file (``CORG.gtf``)
- hisat index files (``CORG.?/ht2``).
- kallisto index file (``CORG.transcripts.Kallisto.indexed``)
- salmon index directory (``CORG.transcripts.Salmon.indexed/``)
- ``COMMANDS.sh`` explaining how each of these files were generated

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

The ``LocalRunExample`` directory contains three `FASTQ <https://en.wikipedia.org/wiki/FASTQ_format>`__ files for CORG containing RNA-seq data. These are examples of local unpaired data, and are each in a directory of their own. The file naming format for these reads is "?\_1.fastq" where the "?" is the number of the sample. GEMmaker finds these files through the glob pattern defined by ``local_samples_path`` in ``nextflow.config``.

The ``RemoteRunExample`` directory contains the file ``SRA_IDs.txt`` which contains a list of names for remote files to be downloaded by GEMmaker from `NCBI's SRA repository <https://www.ncbi.nlm.nih.gov/sra>`__. In this case, there is only one run ID.

Explanation of the Outputs
~~~~~~~~~~~~~~~~~~~~~~~~~~

Once executed, the example should create a directory called ``output`` with several subdirectories. Four of these directories correspond to each sample (3 local, 1 remote), and each of these contains the files generated for that sample. The other directories are the ``reports`` directory and the ``GEMs`` directory.

In each sample directory you will find the following files:

- ``fastq``: The fastq reads file for the experiment.
- ``fastqc``: 6 or 12 files (depending on paired or unpaired data) from fastqc. FastQC is configured to check files before and after trimmomatic.
- ``bam``: Binary alignment file.
- ``ga``: Expression level transcript abundance.
- ``fpkm``: Two-column version of the ``ga`` file with only gene and FPKM value.
- ``tpm``: Two-column version of the ``ga`` file with only gene and TPM value.

The ``reports`` directory will contain a ``multiqc_report.html`` file that provides several statistics about the run.

.. figure:: /images/MultiQC_Report.png
  :alt: MultiQC_Report

Figure 1: Image of the start of the report for the example run when run with Hisat2.

The ``GEMs`` directory contains the final gene-expression matrices (GEMs) in raw, TPM and FPKM form. These GEMs can be used as input to other analyses such as WGCNA and KINC. They can also be visualized as heatmaps -- the heatmap below consists of the FPKM values (divided by 1000) from the local examples. We can see that ``gene_Zeta`` remained constant across all three samples, ``gene_Beta`` decreased, and ``gene_Alpha`` increased.

.. figure:: /images/heatmap.png
  :alt: heatmap

Figure 2: Heatmap of FPKM values from local samples.
