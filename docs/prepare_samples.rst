Step 2: Prepare Sample Data
---------------------------
GEMmaker is capable of processing both locally stored RNA-seq files and automatically downloading samples stored in the `NCBI SRA <https://www.ncbi.nlm.nih.gov/sra>`__ database.  You can provide both types of files to be included in a single run of GEMmaker, or use only local or only remote files.

Using Samples From NCBI SRA
```````````````````````````
GEMmaker supports automatic download and processing of samples from the `NCBI SRA repository <https://www.ncbi.nlm.nih.gov/sra>`__.  To use samples from the SRA, you must first find the list of NCBI SRA Run IDs of the samples you want to process. The run IDs typically start with an **SRR**, **ERR**, or **DRR** prefix.  Do not confuse these with the Experiment IDs which typically start with SRX, ERX or DRX.  The run IDs must be placed, one per line, in a file.

Example of a remote ID File:

.. code:: bash

  SRR1058270
  SRR1058271
  SRR1058272
  SRR1058273
  SRR1058274
  SRR1058275
  SRR1058276
  SRR1058277



Using Samples Stored Locally
````````````````````````````
By default, GEMmaker expects that FASTQ files are uncompressed (not GZ compressed). They can be stored in any directory on the local filesystem.

Paired FASTQ files
''''''''''''''''''
By default, paired files must have a ``_1.fastq`` and a ``_2.fastq`` suffix at the end of the filename.  GEMmaker uses the ``_1`` and ``_2`` designation to differentiate and match paired files.

Non-Paired FASTQ files
''''''''''''''''''''''
By default, if your data is non-paired, GEMmaker expects all files to have a ``_1.fastq`` suffix at the end of the filename.
