.. _examples:

Test GEMmaker
-------------

Quick Test: See How it Works
''''''''''''''''''''''''''''

GEMmaker provides example data to quickly show how it works. This data consists of a small set of local files (contained with GEMmaker) and a remote sample from the `NCBI's SRA repository <https://www.ncbi.nlm.nih.gov/sra>`__. These samples are small to demonstrate usage for a mixed set of local and remote files.  This example assumes you have `Singularity <https://sylabs.io/>`__ installed.

.. note::

    For the examples on this page, Singularity will be used.  Singularity will automatically retrieve the GEMmaker Docker images and by default will store them in the ``work`` folder that Nextflow creates. However, Nextflow may warn that a cache directory is not set. If you intend to run GEMmaker multipe times, you may wish to designate a permanent cache directory by seting the ``NXF_SINGULARITY_CACHEDIR`` prior to running GEMmaker. You can learn more at the `nf-core tools page <https://nf-co.re/tools/#singularity-cache-directory>`_

You can run the example by executing the following command within the GEMmaker directory:

.. code:: bash

   nextflow run systemsgenetics/gemmaker -profile test,<docker/singularity/podman/shifter/charliecloud/conda/institute>


Replace the text `<docker/singularity/podman/shifter/charliecloud/conda/institute>` with the execution profile of your choice. For example to test GEMmaker in a Singularity image the command would be.

.. code:: bash

   nextflow run systemsgenetics/gemmaker -profile test,singularity

Results are stored in the ``results`` directory.

**About the Demo Test Data**

The demo data provided by GEMmaker belongs to the imaginary CORG organism. For the local example, we use a set of 3 artificially made RNA-seq runs. The fictitious CORG organism has a very small "genome" of only 2,336 nucleotides, 3 "chromosomes" and 6 "genes". The 6 genes are named ``gene_Alpha``, ``gene_Beta``, ``gene_Zeta``, ``gene_Gamma``, ``gene_Delta``, ``gene_Epsilon``.

For the remote data file, GEMmaker automatically downloads a very small RNA-seq file from NCBI. This dataset is from an uncharacterized bacteria, but luckily, CORG shares 3 of the genes with this bacteria so we can use CORG's reference file. This remote sample was selected becasue it is an unusually small file, making it  ideal for the example dataset.

Functional Test With Real Data
''''''''''''''''''''''''''''''
The demo data described above is provided to give a new potential user a very quick example for how to use GEMmaker, but the samples and results are meaningless.  If you would like to try using GEMmaker with a real dataset we have included command-line examples in the step-by-step instructions that follow for a set of eight RNA-seq samples available on NCBI aligned to the `TAIR 10 Arabidopsis  genome reference from Ensembl Plants <https://plants.ensembl.org/Arabidopsis_thaliana/Info/Index>`_. The steps will walk you through downloading the reference genome (or transcriptome reference), indexing it with the tool of your preference, and running GEMmaker.

Full Test With Large Data
'''''''''''''''''''''''''
If you would like to test GEMmaker on a large data then the NCBI SRA Project `PRJNA301554 <https://www.ncbi.nlm.nih.gov/bioproject/PRJNA301554/>`_ offers a good choice.  To test this data you can use the following:

- Reference data for Hisat2:

  - *Oryza sativa* `MSU v7.0 genome reference FASTA file <http://rice.uga.edu/pub/data/Eukaryotic_Projects/o_sativa/annotation_dbs/pseudomolecules/version_7.0/all.dir/all.con>`_.
  - *Oryza sativa* `MSU v7.0 annotations in GFF3 <http://rice.uga.edu/pub/data/Eukaryotic_Projects/o_sativa/annotation_dbs/pseudomolecules/version_7.0/all.dir/all.gff3>`_.

- Reference data for Kallisto:

  - *Oryza sativa* `MSU v7.0 transcriptome (cDNA) FASTA file <http://rice.uga.edu/pub/data/Eukaryotic_Projects/o_sativa/annotation_dbs/pseudomolecules/version_7.0/all.dir/all.cdna>`_.
- SRA List of :download:`475 Run IDs <./SRA_IDs.txt>` from the the NCBI SRA Project PRJNA301554.

Follow each step in the following sections and adjust the command-line instructions appropriately for this data.
