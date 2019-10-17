.. _configuration:

The Configuration File
----------------------

GEMmaker uses the ``nextflow.config`` file which allows you to be customize where input data is housed, which tools to execute, which data to publish, and to specify resource limitations for your compute environment. The config file has three main sections:

- ``project``:  Parameters providing background information about the GEMmaker run to be performed.
- ``params``: Parameters for input files, output files, and software.
- ``profiles``: Example profiles for running on different environments such as an HPC system.

The following provides detailed information on each parameter in ``nextflow.config``. Refer to the `Nextflow documentation <https://www.nextflow.io/docs/latest/config.html>`__ for more information on the language of the config file.


Project
~~~~~~~

name
====

A human readable name for the analysis that you will perform. It is best to keep this short and brief.

machine_name
============

A machine readable name for the analysis that you will perform. It is mean to be short. It should only have alphanumeric characters (0-9, a-z, A-Z) and underscores. This parameter will be used to name the final GEM file.

description
===========

In order to help others understand the purpose for the GEMmaker run, you should include a brief description providing enough details to help your collegues who may look at your GEMmaker run in the future.



Input
~~~~~
The input section of the configuration file is for you to specify two types of input data: the genomic reference and the set of input RNA-seq data files.  The type and format of these data will differ depending if you want to use Hiast2, Kallisto or Salmon.   By default, GEMmaker has provided a set of directories to make it easy to organize the input files.  You can find the following directories and files in GEMmaker:

.. code::

  input/
  input/references/
  input/SRA_IDs.txt
  input/samples2skip.txt

You can place your genomic reference index files and the RNA-Seq data files in these directories.  

**Input RNA-Seq Data**

If you are using remote RNA-seq data from NCBI SRA you can edit the ``input/SRA_IDs.txt`` file and add the SRR accession numbers (i.e. run IDs) to this file (one per line).  Alternatively, if you are using local FASTQ files you can place those files directly into the ``input`` directory.  If you have a combination of remote and local files you can do both.  

**Input Genome Reference Data**

All of your genome references index files should go into the ``input/references`` directory.  The demo data that comes with GEMmaker is 
organized into subdirectories and files such as :

.. code::

  CORG.genome.Hisat2.indexed/
  CORG.transcripts.Salmon.indexed/
  CORG.transcripts.Kallisto.indexed

You can follow a similar organization or simply place the index files into the ``input_references`` directory.  Be sure to set the appropriate ``hisat2.index_dir``, ``salmon.index_dir``, or ``kallisto.index_file``.   If you simply place the index files for hisat2 and salmon into the ``input_referecnes`` directory then you can just set the ``hisat2.index_dir`` or ``salmon.index_dir`` to a ``.`` (period).  This tells GEMmaker they are not in a subdirectory.  For the ``kallisto.index_file`` you must always indicate the file name, if it's directly in the ``input_references`` directory or the subdirectory and filename if you have it in a subdirectory.

The following sections describe in more detail the arguments of the Input section of the configuration file.

.. warning::

  GEMmaker provides sample data in the default ``input`` directory to make it easier for someone to test. Before you begin with your own project data, remember to remove the sample data from the ``input`` and ``input/references`` directories.

reference_name
==============
The unique name for the genome reference assembly. It must not contain spaces or special characters, only alphanumeric characters (0-9, a-z, A-Z) and underscores. This name will be used when creating intermediate files that you may want to keep, such as BAM files. 

If Hisat2 is being used as the quanitifaction tool, this name must match the reference name used while running Hisat2 for building indexes in the previous step.

reference_dir
=============
The path to the directory where the genome reference files are housed.  Each quantification tool (Hisat2, Salmon or Kallisto) requires different files.  If no preceeding `/` is used the path is expected to be in the GEMmaker directory. The default is set to the ``input/references`` directory of GEMmaker.

input_data_dir
==============
The path to the directory where any locally stored FASTQ files are housed.  If no preceeding ``/`` is used the path is expected to be in the GEMmaker directory. The default is set to the ``input`` directory of GEMmaker.

remote_sample_list
==================
The path to the file containing a list of SRA Run IDs. These runs will be downloaded from NCBI. This must be a text file with one SRR/DRR/ERR ID per line. No blank lines are allowed. For example:

.. code:: bash

  SRR360147
  SRR493289
  SRR1696865
  SRR2086505
  SRR2086497
  SRR1184187
  SRR1184188

If no remote files are to be downloaded, set this parameter to ``"none"``.  This file must be found in the directory specified by the ``params.input.input_data_dir``.


local_sample_files
==================

The `GLOB <https://en.wikipedia.org/wiki/Glob_(programming)>`__ pattern, that identifies locally stored FASTQ files in the directory specified by the ``input.input_data_dir`` parameter. The default GLOB pattern can find paired or non-paired data that have a ``_1.fastq`` and a ``_2.fastq`` file suffix using the GLOB pattern:

.. code::

  "*_{1,2}.fastq"


hisat2
======

If you want to use the Hisat2 pipeline for alignment and quantification of reads, set ``enable`` to ``true``.   If Hisat2 is enabled, the trimmomatic, samtools and stringtie processes will be enabled as well.

The ``index_dir`` should be the location where the Hisat2 `.ht2` files are located.  This folder must be inside the directory specified by the ``input.reference_dir`` setting.

The ``gtf_file`` parameter should be the name of the GTF file. The file must be inside the directory specified by the ``input.reference_dir`` setting.

Default values:

.. code::

  hisat2 {
      enable = false
      index_dir = "CORG.genome.Hisat2.indexed"
      gtf_file = "CORG.transcripts.gtf"
  }


salmon
======

If you want to use Salmon for quantification of reads, set ``enable`` to ``true``.

The ``index_dir`` should be the name of the directory where Salmon index files are found. These indexes should have been built with from the reference transcript FASTA file using the ``salmon index`` program. The directory must be inside the directory specified by the ``input.reference_dir`` setting.

.. code:: bash

  salmon {
    enable = false
    index_dir = "CORG.transcripts.Salmon.indexed"
  }

kallisto
========

If you want to use Kallisto for quantification of reads, set ``enable`` to ``true``.

The ``index_file`` should be the name of the index file.  This index file should have been built with from the reference genome using the ``kallisto index`` program.  The directory must be inside the directory specified by the ``input.reference_dir`` setting.

.. code:: bash

  kallisto {
    enable = true
    index_file = "CORG.transcripts.Kallisto.indexed"
  }

.. warning::

  You can enable only a Hisat, Kallisto or Salmon but not more than one.



Output
~~~~~~
By default, GEMmaker will store all results in an ``output`` directory that can be found in the GEMmaker directory after GEMmaker runs. This will include several sub directories:

  - sample directories: each sample will have a unique directory with all relevant intermediate files, metadata and log files.
  - ``GEMs``:  will conain the Gene Expression Matricies (GEMs)
  - ``reports``:  will contain MulitQC quality contorl reports.

The output section of the configuration file therefore provides control for where results are saved and which intermediate files should be kept.

.. note::

  The average user will NOT need to change any of the default output parameters.

The following settings and their defaults are :

.. code::

  output {

    // Universal output parameters
    dir = "output"
    sample_dir = { "${params.output.dir}/${sample_id}" }
    publish_mode = "link"
    publish_sra = false
    publish_downloaded_fastq = false
    publish_tpm = true
    publish_raw = true
    multiqc = true
    create_gem = true

    // Salmon and Kallisto specific parameters
    publish_gene_abundance = false

    // Hisat2 specific parameters
    publish_stringtie_gtf_and_ga = false
    publish_trimmed_fastq = false
    publish_bam = false
    publish_sam = false
    publish_fpkm = true
  }

dir
===

All results and reports generated by nextflow are stored in a single output directory.  By default this is set to the ``output`` directory inside of GEMmaker.


sample_dir
==========

Results generated by this workflow are stored in sub directories that are named by their sample ID. If the FASTQ file is not associated with a sample ID (for example, with local files), then the "sample ID" is simply the base name of the FASTQ file.

The default is to have one directory for each sample. However, if you have a large amount of samples (1000s or more), it may be problematic to have so many sample directories in one place. To deal with this issue you can use a pattern that organizes the results into a multi-level directory tree. For example:

.. code:: bash

  sample_dir = { "${params.output.dir}/${sample_id[0..2]}/${sample_id[3..4]}/${sample_id.drop(5)}/${sample_id}" }

This pattern will organize sample directories into three levels of subdirectories. For example, the output of the sample ``SRX0123456`` would be put in the directory ``SRX/12/34/56/SRX123456/``. You can modify the above patterns for your needs.

.. note::

  The pattern shown for the ``sample_dir`` is not a GLOB pattern. It is understood negatively by Nextflow.  The brackets in this example denote a `closure`, a language construct in Nextflow which allows you to create more dynamic expressions using variables and even other configuration parameters. In this case, ``sample_id`` is a variable that will be defined, when GEMMaker runs, for each sample.

publish_mode
============

This controls how intermeidate files are saved.  Options are the standard Nextflow options:

- ``"link"``: Recommended, creates a hardlink for each published file.
- ``"rellink"``: Use when hardlink is not possible.
- ``"symlink"``: Use when hardlink is not possible (currently not compatible with iRODS).
- ``"copy"``: Not recommended, copies each published file to ``publshDir`` after it is created in the pipeline. This option may slow the pipeline significantly.

Intermediate Files
==================

The remaining options in the output parameter determine which intermediate and final output files should be published. By default, all intermediate files are set to false, while final output files are set to true. The following table is a summary of each file:

.. list-table:: Title
   :widths: 25 25 25 50
   :header-rows: 1

   * - Parameter
     - Default
     - Used by
     - Brief Description
   * - publish_sra
     - false
     - Hisat2, Salmon, Kallisto
     - Downloaded Sequence Read Archive (sra) file from NCBI (not human readable)
   * - publish_downloaded_fastq
     - false
     - Hisat2, Salmon, Kallisto
     - Extracted sra file in fastq format (human readable)
   * - publish_tpm
     - true
     - Hisat2, Salmon, Kallisto
     - Transcripts Per Kilobase Million, Final Output Count file option `Extended Descripion <https://www.rna-seqblog.com/rpkm-fpkm-and-tpm-clearly-explained/>`__
   * - publish_raw
     - true
     - Hisat2, Salmon, Kallisto
     - Final Output Count file option, the raw count of each gene. Compare to FPKM and TPM
   * - multiqc
     - true
     - Hisat2, Salmon, Kallisto
     - A final report that is generated that tells you about the GEMmaker run
   * - create_gem
     - true
     - Hisat2, Salmon, Kallisto
     - Combines Final Count Files (FPKM, TPM, raw) into their respective GEM
   * - publish_gene_abundance
     - false
     - Salmon, Kallisto
     - File Generated by Kallisto or Salmon before it is cleaned into Final Count Files
   * - publish_stringtie_gtf_and_ga
     - false
     - Hisat2
     - File Generated by Hisat2 before it is cleaned into Final Count Files
   * - publish_trimmed_fastq
     - false
     - Hisat2
     - Fastq files after they have been trimmed
   * - publish_bam
     - false
     - Hisat2
     - binary alignment file (not human readable) of genes aligned to reference genome
   * - publish_sam
     - false
     - Hisat2
     - alignment file (human readable) of genes aligned to reference genome
   * - publish_fpkm
     - true
     - Hisat2
     - Fragments Per Kilobase Million, Final Output Count file option `Extended Descripion <https://www.rna-seqblog.com/rpkm-fpkm-and-tpm-clearly-explained/>`__



Execution
~~~~~~~~~

queue_size
==========

The maximum number of processes to execute at once.  This is purposely set as a default of 4 to prevent GEMmaker from overrunning a local machine. By default only 4 jobs can execute at a time.  Increase this value appropriate for your local or HPC system resources.

Default:

.. code:: bash

  queue_size = 4


Software
~~~~~~~~
This section is meant to provide customized settings for a software tool. Currently the only tool that requires this is Trimmomatic and Trimmomatic is only used if Hisat2 is enabled.

trimmomatic
===========

The trimmomatic settings and defaults are as follows.

Default:

.. code:: bash

  trimmomatic {
    clip_path = "${baseDir}/files/fasta_adapter.txt"
    MINLEN = "0.7"
    quality = ""
    SLIDINGWINDOW = "4:15"
    LEADING = "3"
    TRAILING = "6"
  }

You should not need to adjust the ``clip_path`` directory unless you have manually installed trimmomatic. If you are using Docker or Singularity with GEMmaker this value show stay as is.  For all others. Please consult the `Trimmomatic documentation <http://www.usadellab.org/cms/?page=trimmomatic>`__ to change these defaults.

Other sections
~~~~~~~~~~~~~~
You will see the following sections present in the configuration file:  ``report``, ``timeline``, ``trace``, ``docker``, ``singularity`` and ``process``.  You should not need to change anything in these sections. To learn more about how they are used, please consult the `Nextflow documentation <https://www.nextflow.io/docs/latest/index.html>`__.

Profiles
~~~~~~~~

The configuration file provides several profiles for running GEMmaker in different computing environments. Each profile defines various settings that override the defaults provided by the rest of the file. The profile that is used by GEMmaker is specified on the command-line at run-time, and they can be combined with each other. For example, to run GEMmaker with the ``pbs`` and ``testing`` profiles enabled:

.. code:: bash

  nextflow run main.nf -profile pbs,testing

You can modify these config files to suit your needs, or even create your own. For more information, refer to the `Nextflow documentation <https://www.nextflow.io/docs/latest/config.html#config-profiles>`__ on config profiles. Here we describe each of the profiles provided by GEMmaker:

docker
======

The ``docker`` profile enables GEMmaker to run processes in Docker containers. This behavior can also be enabled by specifying ``-with-docker`` on the command-line.

k8s
===

The ``k8s`` profile provides basic execution settings for running GEMmaker on a Kubernetes cluster.

modules_kamiak
==============

In lieu of using Docker or Singularity, software dependencies can be provided by environment modules (or a compatible equivalent such as Lmod). Module names tend to vary from system to system. The ``modules_kamiak`` profile is specific to the Washington State University Kamiak cluster. You will likely need to create your own profile that uses the correct module names for your cluster.

modules_palmetto
================

In lieu of using Docker or Singularity, software dependencies can be provided by Environment Modules (or a compatible equivalent such as Lmod). Module names tend to vary from system to system. The ``modules_kamiak`` profile is specific to the Clemson University Palmetto cluster, but you will likely need to create your own profile that uses the correct module names for your cluster.

pbs
===

The ``pbs`` profile provides basic execution settings for running GEMmaker on an HPC system that uses the PBS scheduler. This profile is optimized for the Palmetto cluster at Clemson University, so it may need to be modified to suit your particular system.

singularity
===========

The ``singularity`` profile enables GEMmaker to run processes in Singularity containers. This behavior can also be enabled by specifying ``-with-singularity`` on the command-line.

slurm
=====

The ``slurm`` profile provides basic execution settings for running GEMmaker on an HPC system using the SLURM scheduler. This profile is optimized for the Kamiak cluster at Washington State University, so it may need to be modified to suit your particular system.

standard
========

The ``standard`` profile uses the local executor, in which processes are simply
launched as normal processes on the local machine. By default the local
executor uses the number of CPU cores to limit how many processes are run
in parallel.

testing
=======

The ``testing`` profile overrides the default ``errorStrategy`` to terminate the entire workflow if any error occurs, rather than ignore failed samples. This profile is useful for debugging issues with the workflow, so that the workflow terminates immediately if any process fails.

travis
======

The ``travis`` profile is used by Travis CI for testing new builds.
