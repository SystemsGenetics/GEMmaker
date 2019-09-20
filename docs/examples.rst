.. _examples:

Running a test of GEMmaker
--------------------------

GEMmaker provides example data to quickly show how it works. This data consists of a small set of local files (contained with GEMmaker) and a remote sample from the `NCBI's SRA repository <https://www.ncbi.nlm.nih.gov/sra>`__. These samples are small to demonstrate usage for a mixed set of local and remote files.  This example assumes you have `Singularity <https://sylabs.io/>`__ installed.

You can run the example by executing the following command within the GEMmaker directory:

.. code:: bash

  nextflow run main.nf -profile standard -with-singularity -with-report -with-timeline

The following describes the meaning of the arguments in the command-line above:

   - `-profile standard`: instruct Nextflow to run this example on your local machine
   - `-with-singularity`: instructs Nextflow to use the singularity
   - `-with-report`: instructs Nextflow to generate an HTML report indicating performance of the workflow.
   - `-with-timeline`:  instructs Nextflow to generate an HTML report showing a timline of when steps were executed.

.. note::
  GEMmaker will automatically download the remote files from NCBI SRA and will use Kallisto for transcript expression quantification.

You should see an output that looks like this:

.. code:: bash

  N E X T F L O W  ~  version 19.07.0
  Launching `main.nf` [nasty_heisenberg] - revision: 5d2056118e

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
  Use Hisat2:                 false
  Hisat2 Index Directory:     /local/Projects/GEMmaker/examples/reference/CORG.transcripts.Hisat2.indexed/
  Hisat2 Index Prefix:        CORG
  Hisat2 GTF File:            /local/Projects/GEMmaker/examples/reference/CORG.gtf

  Use Kallisto:               true
  Kallisto Index File:        /local/Projects/GEMmaker/examples/reference/CORG.transcripts.Kallisto.indexed

  Use Salmon:                 false
  Salmon Index File:          /local/Projects/GEMmaker/examples/reference/CORG.transcripts.Salmon.indexed


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
  Queue size:                 4


  Software Parameters:
  --------------------
  Trimmomatic clip path:      /local/Projects/GEMmaker/files/fasta_adapter.txt
  Trimmomatic minimum ratio:  0.7

  executor >  local (37)
  [59/547779] process > retrieve_sra_metadata (1)               [100%] 1 of 1 ✔
  executor >  local (37)
  [59/547779] process > retrieve_sra_metadata (1)               [100%] 1 of 1 ✔
  executor >  local (37)
  [59/547779] process > retrieve_sra_metadata (1)               [100%] 1 of 1 ✔
  executor >  local (37)
  [59/547779] process > retrieve_sra_metadata (1)               [100%] 1 of 1 ✔
  executor >  local (38)
  [59/547779] process > retrieve_sra_metadata (1)               [100%] 1 of 1 ✔
  executor >  local (38)
  [59/547779] process > retrieve_sra_metadata (1)               [100%] 1 of 1 ✔
  executor >  local (38)
  [59/547779] process > retrieve_sra_metadata (1)               [100%] 1 of 1 ✔
  [e7/d57ab2] process > write_stage_files (SRX218012)           [100%] 4 of 4 ✔
  [28/86e2ec] process > start_first_batch                       [100%] 1 of 1 ✔
  [71/38c35e] process > read_sample_file (SRX218012.sample.csv) [100%] 4 of 4 ✔
  [65/e3e664] process > next_sample (4)                         [100%] 4 of 4 ✔
  [8d/6b2215] process > download_runs (SRX218012)               [100%] 1 of 1 ✔
  [21/4c8efd] process > fastq_dump (SRX218012)                  [100%] 1 of 1 ✔
  [a5/35232d] process > fastq_merge (SRX218012)                 [100%] 1 of 1 ✔
  [ab/012e41] process > fastqc_1 (SRX218012)                    [100%] 4 of 4 ✔
  [6b/cf66e1] process > kallisto (SRX218012)                    [100%] 4 of 4 ✔
  [dd/d4bd75] process > kallisto_tpm (SRX218012)                [100%] 4 of 4 ✔
  [-        ] process > salmon                                  -
  [-        ] process > salmon_tpm                              -
  [-        ] process > trimmomatic                             -
  [-        ] process > fastqc_2                                -
  [-        ] process > hisat2                                  -
  [-        ] process > samtools_sort                           -
  [-        ] process > samtools_index                          -
  [-        ] process > stringtie                               -
  [-        ] process > hisat2_fpkm_tpm                         -
  [db/183534] process > multiqc                                 [100%] 1 of 1 ✔
  [2d/2ebe04] process > create_gem                              [100%] 1 of 1 ✔
  [b7/3b52e3] process > clean_sra (SRX218012)                   [100%] 1 of 1 ✔
  [83/45cd5a] process > clean_downloaded_fastq (SRX218012)      [100%] 1 of 1 ✔
  [a4/29eded] process > clean_merged_fastq (SRX218012)          [100%] 1 of 1 ✔
  [-        ] process > clean_trimmed_fastq                     -
  [-        ] process > clean_sam                               -
  [-        ] process > clean_bam                               -
  [78/755446] process > clean_kallisto_ga (SRX218012)           [100%] 4 of 4 ✔
  [-        ] process > clean_salmon_ga                         -
  [-        ] process > clean_stringtie_ga                      -



Results are stored in the ``output`` directory. See the **Explanation of Outputs** section below for information about these files.

.. code:: bash

  output/
    1/
    2/
    3/
    GEMs/
    reports/
    SRX218012/

About the Example Data
~~~~~~~~~~~~~~~~~~~~~~

The example data provided here belongs to the imaginary organism "Cool Organism" (CORG). For the local example, we use a set of 3 artificially made RNA-seq runs. The fictitious CORG organism has a very small "genome" of only 2,336 nucleotides, 3 "chromosomes" and 6 "genes". The 6 genes are named ``gene_Alpha``, ``gene_Beta``, ``gene_Zeta``, ``gene_Gamma``, ``gene_Delta``, ``gene_Epsilon``.

For the remote data file, GEMmaker automatically downloads a very small RNA-seq file from NCBI. This dataset is from an uncharacterized bacteria, but luckily, CORG shares 3 of the genes with this bacteria so we can use CORG's reference file. This remote sample was selected becasue it is an unusually small file, making it  ideal four the example dataset.

Using Salmon or Hisat2
~~~~~~~~~~~~~~~~~~~~~~

By default, GEMmaker uses Kallisto for transcript expression level quantification. If you would like to use Salmon or Hisat2 instead, you must edit ``nextflow.config`` and enable Salmon or Kallisto. In the GEMmaker directory, edit the ``nextflow.config`` file using your favorite text editor. Here we use `vim <https://www.vim.org/>`__ on the command line:

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

Note the ``index_file`` or ``index_dir`` parameters. Each tool uses its own format for indexing the genomic reference. This helps improve speed.  For the example data these indexes already exist.  Therefore, you can save the file and run the worklow using the same command-line as shown previously.


Explanation of the Inputs
~~~~~~~~~~~~~~~~~~~~~~~~~

The inputs for the example run are in the ``input`` directory, and consist of the ``input/references`` directory and all other files for the local and remote samples.

Genome Reference Files
**********************

Hisat2, Kallisto and Salmon use a genome sequence, or reference. Each tool uses its own set of indexes and files. You can find all necessary files for the example CORG genome in the ``input/references/`` directory.

This directory contains:

- ``CORG.fna``: the reference genome file.
- ``CORG.gtf``:  the `GTF <https://uswest.ensembl.org/info/website/upload/gff.html>`__ file listing the gene annotations.
- ``CORG.*.ht2``: A series of Hisat2 index files with the suffix ``ht2``
- ``CORG.transcripts.Kallisto.indexed``: the Kallisto index file.
- ``CORG.transcripts.Salmon.indexed/``: A directory containing Salmon index files.
- ``COMMANDS.sh`` A BASH script with exact commands for creating the indexes.

These are the files needed to run Hisat2, Kallisto, and Salmon on the CORG data.

Sample Data
***********

GEMmaker expects to find all sample data in the ``input`` directory.  Here are three  `FASTQ <https://en.wikipedia.org/wiki/FASTQ_format>`__ files for the local CORG samples. These are examples of local unpaired data. The file naming format for these reads is "?\_1.fastq" where the "?" is the number of the sample. GEMmaker finds these files through the glob pattern defined by ``local_samples_path`` in the ``nextflow.config`` file.

Samples that are found on the NCBI SRA are found in the ``SRA_IDs.txt`` file. This file should contain a list of SRA RUN IDs (i.e. begin with SRR, ERR or DRR) for each sample to be downloaded by GEMmaker from `NCBI's SRA repository <https://www.ncbi.nlm.nih.gov/sra>`__. For this example, there is only one run ID.

Explanation of the Outputs
~~~~~~~~~~~~~~~~~~~~~~~~~~

Once executed, the example should create a directory named ``output``. It will contain a directory for results from each sample: local samples are named `1`, `2`, `3` and the remote sample is named 'SRX218012'.  By default, to save storage space, GEMmaker will only place log files or analysis reports for each sample.  Although, you can choose to have GEMmaker include downlaoded FASTQ files (for remote samples), trimmed FASTQ and BAM (if Hisat2 is used), or abundance files (if Salmon or Kallisto are used).

Additionally, a ``reports`` directory containing the `MultiQC <https://multiqc.info/>`__ summary of performance for the bioinformatics tool:


.. figure:: /images/MultiQC_Report.png
  :alt: MultiQC_Report

Figure 1: Image of the start of the report for the example run when run with Hisat2.

The ``GEMs`` directory contains the final gene-expression matrices (GEMs) in raw, TPM and FPKM form. These GEMs can be used as input to other analyses such as WGCNA and KINC. They can also be visualized as heatmaps -- the heatmap below consists of the FPKM values (divided by 1000) from the local examples. We can see that ``gene_Zeta`` remained constant across all three samples, ``gene_Beta`` decreased, and ``gene_Alpha`` increased.

.. figure:: /images/heatmap.png
  :alt: heatmap

Figure 2: Heatmap of FPKM values from local samples.
