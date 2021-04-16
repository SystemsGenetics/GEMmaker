.. _execution:

Step 3: How to Run GEMmaker
---------------------------

Configuration
'''''''''''''
The instructions here provide details for running GEMmaker using Singularity.
For most instances you probably won't need to make customizations to the workflow configuration. However, should you need to, GEMmaker is a `nf-core <https://nf-co.re/>`_ compatible workflow.  Therefore, it follows the general approach for workflow configuration which is described at the `nf-core Pipeline Configuration page <https://nf-co.re/usage/configuration>`_.  Please see those instructions for the various platforms and settings you can configure.  You will find GEMmaker specific configuration options here.

Quick Start
'''''''''''
The following examples for running GEMmaker use the Arabidopsis thaliana reference genome available from `Ensembl Plants <https://plants.ensembl.org/Arabidopsis_thaliana/Info/Index>`_. The reference sequences in the examples are prepared the same as described in Step 2.  As an example, we will indicate 3 SRA files for automatic retrieval and processing by listing them in a file named ``SRAs.txt``:

.. code:: bash

    SRR1058270
    SRR1058271
    SRR1058272

Kallisto
........
To run Kallisto you need to specify:

- The path to the genome reference indexed file
- A file containing a set of SRA run IDs you want to download or the path were FASTQ files are stored on the local system.

For example:

.. code:: bash

  nextflow run systemsgenetics/gemmaker -profile singularity \
    --pipeline kallisto \
    --kallisto_index_path Arabidopsis_thaliana.TAIR10.kallisto.indexed \
    --sras SRAs.txt

Salmon
......
To run Salmon you need to specify:

 - The path to the directory containing the genome reference index files.
 - A file containing a set of SRA run IDs you want to download or the path were FASTQ files are stored on the local system.

 For example:

.. code:: bash

  nextflow run systemsgenetics/gemmaker -profile singularity \
    --pipeline salmon \
    --salmon_index_path Arabidopsis_thaliana.TAIR10.salmon.indexed \
    --sras SRAs.txt

Hisat2
......
To run Hiast2 you need to specify:

- The path to directory containing the Hisat2 genome reference indexed files
- The base name of the whole genome. All Hisat2 index files use this base name. For this example, the base name used is  ``Arabidopsis_thaliana.TAIR10``.
- The GTF file containing the gene annotations.
- A file containing a set of SRA run IDs you want to download or the path were FASTQ files are stored on the local system.

For example:

.. code:: bash

  nextflow run systemsgenetics/gemmaker -profile singularity \
    --pipeline hisat2 \
    --sras SRAs.txt \
    --hisat2_base_name Arabidopsis_thaliana.TAIR10 \
    --hisat2_index_dir hisat2_indexes \
    --hisat2_gtf_file Arabidopsis_thaliana.TAIR10.gtf

Additionally, you can control the Trimmomatic trimming step by adding any of the following parameters:

- ``--trimmomatic_clip_file``: the location for a custom file of sequences to clip. GEMmaker provides a default version so you only need to set this if you have custom sequences.
- ``--trimmomatic_MINLEN``: corresponds to the ``MINLEN`` argument of Trimmomatic. Defaults to 0.7.
- ``--trimmomatic_SLIDINGWINDOW``: corresponds to the ``SLIDINGWINDOW`` argument of Trimmomatic. Defaults to "4:15"
- ``--trimmomatic_LEADING``: corresponds to the ``LEADING`` argument of Trimmomatic. Defults to 3.
- ``--trimmomatic_TRAILING``: correponds to teh ``TRAILING`` argument of Trimmomatic. Defaults to 6.

Intermediate Files
''''''''''''''''''
GEMmaker was designed to limit the storage requirements in order to allow for processing of large numbers of FASTQ files without overrunning storage requirement.  By default it will remove all large intermediate files to keep space usage to a minimum. However, you can indicate what intermediate files you would like to keep by providing any of the following arguments and setting them to ``true``.  For example, to keep the downloaded SRA files the ``keep_sra``` argument would be provided and set to true:

.. code:: bash

  nextflow run systemsgenetics/gemmaker -profile singularity \
    --pipeline salmon \
    --salmon_index_path Arabidopsis_thaliana.TAIR10.salmon.indexed \
    --sras SRAs.txt \
    --keep_sra true

The following is a listing of all arguments that can control which intermediate files are kept.

SRA Files
.........
The following arguments can be used if the ``--sras`` option is used.

- ``--keep_sra``: Set to true to keep all downloaded SRA files .
- ``--keep_retrieved_fastq``: Set to true to keep the FASTQ files that are derived from downloaded SRA files.

Kallisto Files
..............
The following arguments can be used if the ``--pipeline kallisto`` option is used.

- ``--kallisto_keep_data``: Set to true to keep the intermediate files created by Kallisto.

Salmon Files
............
The following arguments can be used if the ``--pipeline salmon`` option is used.

- ``--kallisto_keep_data``: Set to true to keep the intermediate files created by Salmon.

Hisat2 Files
............
The following arguments can be used if the ``--pipeline hisat2`` option is used.

- ``--hisat2_keep_data``: Set to true to keep the stringtie output.
- ``--hisat2_keep_sam``: Set to true to keep the SAM files created by Hisat2.
- ``--hisat2_keep_bam``: Set to true to keep the BAM files created by Hisat2.
- ``--trimmomatic_keep_trimmed_fastq``: Set to true to keep the trimmed FASTQ files after trimmomatic is run.
