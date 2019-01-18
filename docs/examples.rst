.. _examples:

Example Usage
-------------

The GEMmaker example consists of a small set of local files (provided in the GEMmaker repo) and a remote file (automatically downloaded from `NCBI's SRA repository <https://www.ncbi.nlm.nih.gov/sra>`__). These samples are extremely small and are only meant to demonstrate usage for a mixed set of local and remote files.

Once GEMmaker and its dependencies have been installed, you can run GEMmaker with the example data. First, copy the example config file to ``nextflow.config``:

.. code:: bash

  cp nextflow.config.example nextflow.config

The ``nextflow.config.example`` file is already configured to run this example. Execute the workflow using the following command:

.. code:: bash

  nextflow run main.nf -profile standard -with-docker

You should see an output that looks like this:

.. code:: bash

  N E X T F L O W  ~  version 18.10.1
  Launching `main.nf` [peaceful_gutenberg] - revision: 137bc7ccff

  ===================================
   G E M M A K E R   P I P E L I N E
  ===================================

  General Information:
  --------------------
    Profile(s):         standard
    Container Engine:   docker


  Input Parameters:
  -----------------
    Remote fastq list path:     /home/{usr}/GEMmaker/examples/RemoteRunExample/SRA_IDs.txt
    Local sample glob:          /home/{usr}/GEMmaker/examples/LocalRunExample/Sample*/\*_{1,2}.fastq
    Reference genome path:      /home/{usr}/GEMmaker/examples/reference/
    Reference genome prefix:    CORG


  Output Parameters:
  ------------------
    Output directory:           /home/{usr}/GEMmaker/output
    Publish downloaded FASTQ:   true
    Publish trimmed FASTQ:      true
    Publish BAM:                true
    Publish FPKM:               true
    Publish TPM:                true


  Execution Parameters:
  ---------------------
    Queue size:                 100
    Number of threads:          1
    Maximum retries:            2
    Error strategy:             ignore


  Software Parameters:
  --------------------
    Trimmomatic clip path:      /home/{usr}/GEMmaker/files/fasta_adapter.txt
    Trimmomatic minimum ratio:  0.7

  [warm up] executor > local
  [98/cba9bc] Submitted process > write_stage_files (1)
  [d1/b1cbfd] Submitted process > write_stage_files (2)
  [fb/ea14ef] Submitted process > write_stage_files (3)
  [a1/822d07] Submitted process > retrieve_sample_metadata
  [10/ff1219] Submitted process > write_stage_files (SRX218012)
  [26/3da51b] Submitted process > start_first_batch
  [6a/9ab954] Submitted process > read_sample_file (1.sample.csv)
  [74/f3836d] Submitted process > read_sample_file (2.sample.csv)
  [1b/d263eb] Submitted process > read_sample_file (3.sample.csv)
  [0a/695151] Submitted process > read_sample_file (SRX218012.sample.csv)
  [13/6fc3ea] Submitted process > fastqc_1 (1)
  [ce/061c10] Submitted process > trimmomatic (1)
  [2a/81acc1] Submitted process > fastqc_1 (2)
  [82/019602] Submitted process > trimmomatic (2)
  [f3/7a17ae] Submitted process > fastqc_1 (3)
  [75/b755a2] Submitted process > fastq_dump (SRX218012)
  [01/4c7c4f] Submitted process > trimmomatic (3)
  [dc/0e3904] Submitted process > fastqc_2 (2)
  [b7/eb221f] Submitted process > hisat2 (2)
  [38/98d3ce] Submitted process > fastqc_2 (1)
  [c0/52f716] Submitted process > hisat2 (1)
  [be/1be15c] Submitted process > SRR_combine (SRX218012)
  [9f/aa1f14] Submitted process > fastqc_2 (3)
  [05/66f559] Submitted process > hisat2 (3)
  [f7/8f2780] Submitted process > fastqc_1 (SRX218012)
  [aa/0d9e57] Submitted process > trimmomatic (SRX218012)
  [94/cedb63] Submitted process > samtools_sort (2)
  [38/454c25] Submitted process > samtools_sort (3)
  [36/f87963] Submitted process > samtools_sort (1)
  [6d/458628] Submitted process > fastqc_2 (SRX218012)
  [cd/05aa92] Submitted process > hisat2 (SRX218012)
  [c3/298c49] Submitted process > samtools_index (2)
  [ed/e775a0] Submitted process > samtools_index (3)
  [7a/57bb71] Submitted process > samtools_index (1)
  [90/62173a] Submitted process > samtools_sort (SRX218012)
  [08/d28f6c] Submitted process > stringtie (3)
  [cf/4e30b3] Submitted process > stringtie (1)
  [c5/f89c37] Submitted process > samtools_index (SRX218012)
  [d7/67724f] Submitted process > stringtie (2)
  [ca/881318] Submitted process > stringtie (SRX218012)
  [fc/5688e8] Submitted process > hisat2_raw (3)
  [30/93eb53] Submitted process > fpkm_or_tpm (3)
  [91/969c3a] Submitted process > hisat2_raw (SRX218012)
  [9b/9c541f] Submitted process > fpkm_or_tpm (SRX218012)
  [49/ddb561] Submitted process > fpkm_or_tpm (1)
  [1b/3dbd3d] Submitted process > hisat2_raw (1)
  [df/c3f00c] Submitted process > fpkm_or_tpm (2)
  [5c/3053f4] Submitted process > hisat2_raw (2)
  [32/df5310] Submitted process > next_sample (1)
  [ea/812195] Submitted process > multiqc
  [9c/d98d23] Submitted process > createGEM

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
