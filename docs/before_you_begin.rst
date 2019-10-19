Before You Begin
----------------
Before you get started there are some important things to know.  To use GEMmaker you must manually prepare the genomic reference data, organize your samples and configure the workflow with settings to deliver the outputs you desire.  Instructions for this are provided in a later section.

.. important ::

  You should only use GEMmaker with one set of samples.  You can use as many samples as you like, but avoid changing the sample list between runs of GEMmaker.  If you have multiple projects each with a different set of samples, you should install an instance of GEMmaker for each project.   Also, if you have used GEMmaker for testing, it is best to not reuse the same GEMmaker installation used for testing with an actual sample set.  Instead create a new installation.

About the Configuration File
````````````````````````````
GEMmaker uses the ``nextflow.config`` file to provide customizations. The configuration file has three primary sections:

- ``project``:  Parameters providing background information about the GEMmaker run to be performed.
- ``params``: Parameters for input files, output files, and software.
- ``profiles``: Example profiles for running on different environments such as an HPC system.

To use GEMmaker you must edit the ``nextflow.config`` file and update the settings appropriate for your data. The file is pre-configured for the example data that comes with GEMmaker. Appropriate values for the settings in the config file will be described in the :doc:`./configuration` section.

Initial Setup
`````````````
Before proceeding, please clear out the example data that comes with GEMmaker. You can do so by executing the following commands from within the GEMmaker directory:

.. code:: bash

  # Remove the example FASTQ files
  rm ./input/*.fastq

  # Remove the example genome reference files
  rm -rf ./inputs/reference/*

  # Clear out the list of  remote NCBI SRA samples to retrieve.
  echo "" > ./inputs/SRA_IDs.txt
