Prepare Sample Data
-------------------
GEMmaker is capable of processing both locally stored RNA-seq files and automatically downloading samples stored in the `NCBI SRA <https://www.ncbi.nlm.nih.gov/sra>`__ database.  You can provide both types of files to be included in a single run of GEMmaker, or use only local or only remote files.

.. warning::

  Before proceeding please remember to clear out the ``input`` directory of example data as described on the :doc:`./usage` page, or you may accidentally include example data with your own data.

Using Samples From NCBI SRA
```````````````````````````
GEMmaker supports automatic download and processing of samples from the `NCBI SRA repository <https://www.ncbi.nlm.nih.gov/sra>`__.  To use samples from the SRA, you must first find the list of NCBI SRA Run IDs of the samples you want to process. The run IDs typically start with an **SRR**, **ERR**, or **DRR** prefix.  Do not confuse these with the Experiment IDs which typically start with SRX, ERX or DRX.  The run IDs must be placed, one per line, in the file ``input/SRA_IDs.txt`` file.

Example of a remote ID File:

.. code:: bash

  SRR360147
  SRR493289
  SRR1696865
  SRR2086505
  SRR2086497
  SRR1184187
  SRR1184188

Not Using Samples From NCBI SRA
```````````````````````````````
If you do not wish to use samples from the NCBI SRA, you must edit the ``input.remote_sample_list``  in the ``nextflow.config`` file and set it to ``"none"``.

Using Samples Stored Locally
````````````````````````````
By default, GEMmaker expects that FASTQ files are uncompressed and stored in the ``input`` directory of GEMmaker.

Paired FASTQ files
''''''''''''''''''
By default, paired files must have a ``_1.fastq`` and a ``_2.fastq`` suffix at the end of the filename.  GEMmaker uses the ``_1`` and ``_2`` designation to differentiate and match paired files.

Non-Paired FASTQ files
''''''''''''''''''''''
By default, if your data is non-paired, GEMmaker expects all files to have a ``_1.fastq`` suffix at the end of the filename.

If your sample naming is different
''''''''''''''''''''''''''''''''''
If you FASTQ files follow some other naming format other than described above for either paired or non-paried data. You can adjust the ``input.local_sample_files`` parameter in the ``nextflow.config`` file.   This parameter requires that you use a `GLOB <https://en.wikipedia.org/wiki/Glob_(programming)>__` pattern to identify files.  For example, the GLOB pattern used by default for GEMmaker to find files with a ``_1.fastq` or ``_2.fastq`` suffix is:

.. code::

  "*_{1,2}.fastq"

If your files are stored elsewhere
''''''''''''''''''''''''''''''''''
As mentioned previously, the default expectation is that all FASTQ files are present in the ``input`` directory of GEMmaker.  If your files are stored elsewhere, you do not have to move them.  Edit the ``nextflow.config`` file and change the ``input.input_data_dir`` parameter to include the full path where the FASTQ files are housed.

If you do not have local files
''''''''''''''''''''''''''''''
If you do not have locally stored files, you must edit the ``input.local_sample_files``  in the ``nextflow.config`` file and set it to ``"none"``.
