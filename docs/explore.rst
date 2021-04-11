.. _examples:

Test GEMmaker
-------------

GEMmaker provides example data to quickly show how it works. This data consists of a small set of local files (contained with GEMmaker) and a remote sample from the `NCBI's SRA repository <https://www.ncbi.nlm.nih.gov/sra>`__. These samples are small to demonstrate usage for a mixed set of local and remote files.  This example assumes you have `Singularity <https://sylabs.io/>`__ installed.

You can run the example by executing the following command within the GEMmaker directory:

.. code:: bash

   nextflow run nf-core/gemmaker -profile test,<docker/singularity/podman/shifter/charliecloud/conda/institute>


Replace the text `<docker/singularity/podman/shifter/charliecloud/conda/institute>` with the execution profile of your choice. For example to test GEMmaker in a Singularity image the command would be.

.. code:: bash

   nextflow run nf-core/gemmaker -profile test,singularity

Results are stored in the ``results`` directory. You can find more information about the results in the **Use GEMmaker section**

About the Demo Test Data
~~~~~~~~~~~~~~~~~~~~~~~~

The demo data provided by GEMmaker belongs to the imaginary CORG organism. For the local example, we use a set of 3 artificially made RNA-seq runs. The fictitious CORG organism has a very small "genome" of only 2,336 nucleotides, 3 "chromosomes" and 6 "genes". The 6 genes are named ``gene_Alpha``, ``gene_Beta``, ``gene_Zeta``, ``gene_Gamma``, ``gene_Delta``, ``gene_Epsilon``.

For the remote data file, GEMmaker automatically downloads a very small RNA-seq file from NCBI. This dataset is from an uncharacterized bacteria, but luckily, CORG shares 3 of the genes with this bacteria so we can use CORG's reference file. This remote sample was selected becasue it is an unusually small file, making it  ideal for the example dataset.
