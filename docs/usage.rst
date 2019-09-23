Use GEMmaker
------------
To use GEMmaker you must prepare the genomic reference data, organize your samples and configure the workflow with settings to deliver the outputs you desire.

Configuration File
``````````````````
GEMmaker uses the ``nextflow.config`` file to provide customizations. The configuration file has three primary sections:

- ``project``:  Parameters providing background information about the GEMmaker run to be performed.
- ``params``: Parameters for input files, output files, and software.
- ``profiles``: Example profiles for running on different environments such as an HPC system.

Settings from this config file will be referred to in the following sections and more details will be provided in the :doc:`./configuration` page.

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


Using GEMmaker
``````````````
The following sections provide step-by-step instructions to prepare, configure and execute GEMmaker.

.. toctree::
   :maxdepth: 1

   prepare_reference
   prepare_samples
   configuration
   execution
