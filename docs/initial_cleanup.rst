Step 1: Prepare GEMmaker
------------------------

Before proceeding, you must clear out the example data that comes with GEMmaker. You can do so by executing the following commands from within the GEMmaker directory:

.. code:: bash

  # Remove the example FASTQ files
  rm ./input/*.fastq

  # Remove the example genome reference files
  rm -rf ./inputs/reference/*

  # Clear out the list of  remote NCBI SRA samples to retrieve.
  echo "" > ./inputs/SRA_IDs.txt

.. note::

  If you have run the demo data as described in the previous section for testing GEMmaker you should use a new installation of GEMmaker for your own data.  You should avoid using the same GEMmaker instance for different sets of samples.
