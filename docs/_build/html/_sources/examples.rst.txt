Running the Examples
--------------------

GEMmaker provides two examples for processing RNA-seq samples that are
either stored on a local computer or housed in `NCBI's SRA
repository <https://www.ncbi.nlm.nih.gov/sra>`__. You can test this
workflow without installing all of the necessary prerequisites by using
the
`GEMmaker-docker <https://github.com/SystemsGenetics/GEMmaker-docker>`__
image.

Analyze Local Data
~~~~~~~~~~~~~~~~~~

This example uses the imaginary organism "Cool Organism" (CORG) and a
data set of 3 artificially made RNA-seq runs for CORG. CORG has a very
small "genome" of only 1,953 nucleotides, 2 "chromosomes" and 3 "genes".
The 3 genes are named "gene\_Alpha", "gene\_Beta" and "gene\_Zeta". This
data set is very small so it can be run on a desktop machine in a short
amount of time.

Reference directory
===================

The reference directory for the is located at:

.. code:: bash

    GEMmaker/examples/LocalRunExample/reference

This directory contains the

- made up reference genome file (**CORG.fna**),
- `GTF <https://uswest.ensembl.org/info/website/upload/gff.html>`__ file (**CORG.gtf**)
- hisat index files (**CORG.?/ht2**).

Data directory
==============

The example data directory is located at:

.. code:: bash

    /GEMmaker/examples/LocalRunExample/Data/

This directory contains 3
`FASTQ <https://en.wikipedia.org/wiki/FASTQ_format>`__ files containing
RNA-seq data. These are examples of unpaired data, and are each in a
directory of their own. The file naming format for these reads is
"?\_1.fastq" where the "?" is the number of the sample. GEM-maker finds
these files through the glob pattern assigned to the
"local\_samples\_path" in the **nextflow.config** file.

Executing the example
=====================

To execute the GEMmaker with the example data you must first rename the
**nextflow.config.example** file as **nextflow.config**:

.. code:: bash

    mv nextflow.config.example nextflow.config

The nextflow.config.example is already configured with the proper
settings to use the local example data just described. To execute the
workflow run the following:

.. code:: bash

    nextflow run main.nf -profile standard

Results
=======

Once executed, the Local Run Example should output 3 directories. These
will be titled Sample\_1, Sample\_2, and Sample\_3. In a real run,
GEMmaker will automatically combine files that have the same experiment
number( [SED]RX0000000 ) but different run numbers ( [SED]RR0000000 ),
so it is possible that the [SED]RX number contains multiple [SED]RR
runs. However, in the the local example, this is not the case.

In each output directory you will find the following files:

- **fastq** The fastq reads file for the experiment.
- **fastqc** 6 or 12 files (depending on paired or unpaired data)
  from fastqc. Fastqc is set up tocheck files before and after trimmomatic.
- **sam** alignment file.
- **bam** binary alignment file.
- **ga** expression level transcript abundance.
- **fpkm** 2 column version of **ga** file with only gene and FPKM value.
- **tpm** 2 column version of the **ga** file with only gene and TPM values.

The output of GEM-maker can be used for several different analysis. The
FPKM files can be combined into an expression matrix and then visualized
using a heatmap. The following heatmap is the Local Example's fpkm
values divided by 1000 in heatmap form. We can see that gene\_Zeta
remained constant across all three samples, gene\_Beta decreased, and
gene\_Alpha increased.

.. figure:: /images/heatmap.png
   :alt: heatmap

   heatmap

Analyze Remote Data
~~~~~~~~~~~~~~~~~~~

This example uses data from an unidentified bacteria. The data set is
located on NCBI, and is unusually small, which makes it useful as an
example. This datasset was selected as an example because of it's small
size which can be run on a desktop computer.

Reference directory
===================

The reference directory containing the genome data for the this example
is located at:

.. code:: bash

    GEMmaker/examples/RemoteRunExample/reference

It contains the reference files like the directory for the local
example, but for the unidentified bacteria rather than CORG.

Input File
==========

Unlike for local data, analysis of remote data requires a list of `NCBI
Sequence Read Archive (SRA) <https://www.ncbi.nlm.nih.gov/sra>`__. Run
IDs. The file should contain one SRA RUN ID per line, with no blank
lines at the bottom. Here is an example file:

.. code:: bash

    GEMmaker/examples/RemoteRunExample/SRA_IDs.txt

The file specifies the SRA RUN ID
`SRR649944 <https://www.ncbi.nlm.nih.gov/sra/SRR649944/>`__

Executing the example
=====================

To run the Remote Run Example, you must change some parameters in the
**nextflow.config** file.

First, change the **params.software\_params.hisat2.path** so that it
points to the proper reference files directory. For this example, it
should be set to:

.. code:: bash

    path = "${PWD}/examples/RemoteRunExample/reference/"

The parameter below this is **params.software\_params.hisat2.prefix**.
This should be changed to the prefix of the reference files. In this
example, the prefix should be set to:

.. code:: bash

    prefix = "GCA_002793175.1_ASM279317v1_genomic"

Additionally, **params.input.remote\_list\_path** should be set to point
to the SRA\_IDs.txt file:

.. code:: bash

    remote_list_path = "${PWD}/examples/RemoteRunExample/SRA_IDs.txt"

The **params.input.local\_samples\_path** paramter should be set to
indicate that there are no local samples:

.. code:: bash

    local_samples_path = "none"

Results
=======

Because there is only one smaple, a single directory named "SRR649944"
will be presetn. This is the NCBI experiment ID that the sample belongs
to. It will contain the same output files describe above for the local
run example. It should be noted that this RNA-seq dataset does not
produce any fpkm values.

Analyzing Both Local and Remote Together
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

GEMmaker is capable of combining both local and remote data in a single
execution if both the **remote\_list\_path** and
**local\_samples\_path** are set in the **nextflow.config** file.
