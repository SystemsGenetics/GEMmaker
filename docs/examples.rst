.. _examples:

Example
-------

The GEMmaker example consists of a small set of local files (provided in the GEMmaker repo) and a remote file (automatically downloaded from `NCBI's SRA repository <https://www.ncbi.nlm.nih.gov/sra>`__). These samples are extremely small and are only meant to demonstrate usage for a mixed set of local and remote files.

Once GEMmaker and its dependencies have been installed, you can run GEMmaker with the example data. The ``nextflow.config`` file is already configured to run this example. Execute the workflow using the following command:

.. code:: bash

  nextflow run main.nf -profile standard -with-singularity

You should see an output that looks like this:

.. code:: bash

  N E X T F L O W  ~  version 19.04.0
  Launching `main.nf` [maniac_roentgen] - revision: db373fdbf9

  ===================================
  G E M M A K E R   P I P E L I N E
  ===================================

  Workflow Information:
  ---------------------
  Project Directory:  /local/Projects/GEMmaker
  Launch Directory:   /local/Projects/GEMmaker
  Work Directory:     /local/Projects/GEMmaker/work
  Config Files:       [/local/Projects/GEMmaker/nextflow.config]
  Container Engine:   singularity
  Profile(s):         standard


  Input Parameters:
  -----------------
  Remote fastq list path:     /local/Projects/GEMmaker/examples/RemoteRunExample/SRA_IDs.txt
  Local sample glob:          /local/Projects/GEMmaker/examples/LocalRunExample/Sample*/*_{1,2}.fastq


  Quantification Tool Input:
  --------------------------
  Use Hisat2:                 true
  Use Kallisto:               false
  Use Salmon:                 false

  Hisat2 Index Directory:     /local/Projects/GEMmaker/examples/reference/CORG.transcripts.Hisat2.indexed/
  Hisat2 Index Prefix:        CORG
  Hisat2 GTF File:            /local/Projects/GEMmaker/examples/reference/CORG.gtf


  Output Parameters:
  ------------------
  Output directory:           output
  Publish SRA:                false
  Publish downloaded FASTQ:   false
  Publish trimmed FASTQ:      false
  Publish BAM:                false
  Publish Gene Abundance:     false
  Publish GTF_GA:             false
  Publish RAW:                false
  Publish FPKM:               true
  Publish TPM:                true
  MultiQC:                    true
  Create GEM:                 false


  Execution Parameters:
  ---------------------
  Queue size:                 100


  Software Parameters:
  --------------------
  Trimmomatic clip path:      /local/Projects/GEMmaker/files/fasta_adapter.txt
  Trimmomatic minimum ratio:  0.7

  executor >  local (69)
  [26/838608] process > write_stage_files      [100%] 4 of 4 ✔
  [08/05fb92] process > retrieve_sra_metadata  [100%] 1 of 1 ✔
  [03/a3e086] process > start_first_batch      [100%] 1 of 1 ✔
  [18/623434] process > read_sample_file       [100%] 4 of 4 ✔
  [50/c71960] process > fastqc_1               [100%] 4 of 4 ✔
  [c6/cba919] process > trimmomatic            [100%] 4 of 4 ✔
  [91/ce105f] process > prefetch               [100%] 1 of 1 ✔
  [c5/d53185] process > fastqc_2               [100%] 4 of 4 ✔
  [e0/3e75af] process > hisat2                 [100%] 4 of 4 ✔
  [37/1ec3af] process > samtools_sort          [100%] 4 of 4 ✔
  [03/fd7e0d] process > samtools_index         [100%] 4 of 4 ✔
  [6b/c29fb9] process > clean_sam              [100%] 4 of 4 ✔
  [a3/b34d17] process > stringtie              [100%] 4 of 4 ✔
  [9e/278762] process > fastq_dump             [100%] 1 of 1 ✔
  [b5/ba3bfb] process > clean_bam              [100%] 4 of 4 ✔
  [da/c708ae] process > hisat2_fpkm_tpm        [100%] 4 of 4 ✔
  [31/803c44] process > clean_trimmed_fastq    [100%] 4 of 4 ✔
  [6c/53f998] process > next_sample            [100%] 4 of 4 ✔
  [87/946d31] process > clean_stringtie_ga     [100%] 4 of 4 ✔
  [52/18bba8] process > fastq_merge            [100%] 1 of 1 ✔
  [78/2e8d98] process > clean_sra              [100%] 1 of 1 ✔
  [24/07979d] process > clean_downloaded_fastq [100%] 1 of 1 ✔
  [02/2a6b86] process > multiqc                [100%] 1 of 1 ✔
  [d8/32b555] process > clean_merged_fastq     [100%] 1 of 1 ✔
  Completed at: 04-Jul-2019 02:29:25
  Duration    : 15.9s
  CPU hours   : (a few seconds)
  Succeeded   : 69





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

The example uses Hisat2 by default. If you would like to use Salmon or Kallisto instead, you must edit ``nextflow.config`` and enable Salmon or Kallisto. In the GEMmaker directory, edit ``nextflow.config`` using your favorite text editor. Here we use `vim <https://www.vim.org/>`__ on the command line:

.. code:: bash

  vim nextflow.config

Then edit the ``params.input`` section in the config file. Set ``enable`` to ``true`` for either Salmon or Kallisto. For example, to use Kallisto:

.. code:: bash

  hisat2 {
    enable = false
    index_dir = "${baseDir}/examples/reference/CORG.transcripts.Hisat2.indexed/"
    index_prefix = "CORG"
    gtf_file = "${baseDir}/examples/reference/CORG.gtf"
  }
  salmon {
    enable = false
    index_dir = "${baseDir}/examples/reference/CORG.transcripts.Salmon.indexed"
  }
  kallisto {
    enable = true
    index_file = "${baseDir}/examples/reference/CORG.transcripts.Kallisto.indexed"
  }

Then save your file and run the worklow:

.. code:: bash

  nextflow run main.nf -profile standard -with-docker

Explanation of the Inputs
~~~~~~~~~~~~~~~~~~~~~~~~~

The inputs for the example run are in the ``examples`` directory, and they consist of the reference directory and two data directories for local and remote samples.

Genome Reference Assembly Files
===============================

Hisat2, Kallisto and Salmon use a genome sequence, or reference. Each tool uses its own set of indexes and files. You can find all necessary files for the example CORG genome here:

.. code:: bash

  GEMmaker/examples/reference/

This directory contains the

- reference genome file (``CORG.fna``),
- `GTF <https://uswest.ensembl.org/info/website/upload/gff.html>`__ file (``CORG.gtf``)
- Hisat2 index files (``CORG.transcripts.Hisat2.indexed/ht2``).
- Kallisto index file (``CORG.transcripts.Kallisto.indexed``)
- Salmon index directory (``CORG.transcripts.Salmon.indexed/``)
- ``COMMANDS.sh`` explaining how each of these files were generated

These are the files needed to run Hisat2, Kallisto, and Salmon on the CORG data.

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
