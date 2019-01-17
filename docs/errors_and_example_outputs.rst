Example Output
==============

When the example data is run using this command:

.. code:: bash
  nextflow run main.nf -profile standard,inDocker

You should see this output on the command line:
.. code:: bash
  N E X T F L O W  ~  version 18.10.1
  Launching `main.nf` [peaceful_gutenberg] - revision: 137bc7ccff

  ===================================
   G E M M A K E R   P I P E L I N E
  ===================================

  General Information:
  --------------------
    Profile(s):         standard,localDocker
    Container Engine:   docker


  Input Parameters:
  -----------------
    Remote fastq list path:     /home/{usr}/GEMmaker/examples/RemoteRunExample/SRA_IDs.txt
    Local sample glob:          /home/{usr}/GEMmaker/examples/LocalRunExample/Sample*/*_{1,2}.fastq
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

Example Files

After the example is run, the output directory should contain these directories:

.. code:: bash
  output/
    1/
    2/
    3/
    GEM/
    reports/
    SRX218012/
