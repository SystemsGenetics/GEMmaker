.. _execution:

Step 3: Run GEMmaker
--------------------

How to Launch GEMmaker
''''''''''''''''''''''
The demonstrate how to to use GEMmaker the `Arabidopsis thaliana` reference genome available from `Ensembl Plants <https://plants.ensembl.org/Arabidopsis_thaliana/Info/Index>`_ was prepared in Step 2.  As an example, we will indicate 3 SRA files for automatic retrieval and processing by listing them in a file named ``SRAs.txt``:

.. code:: bash

    SRR1058270
    SRR1058271
    SRR1058272

If you followed the example in the previous step you should have the reference genome already indexed.

Use Kallisto
............
To run Kallisto you need to specify:

- The path to the genome reference indexed file
- A file containing a set of SRA run IDs you want to download or the path were FASTQ files are stored on the local system.

For example:

.. code:: bash

  nextflow run systemsgenetics/gemmaker -profile singularity \
    --pipeline kallisto \
    --kallisto_index_path Arabidopsis_thaliana.TAIR10.kallisto.indexed \
    --sras SRAs.txt

Use Salmon
..........
To run Salmon you need to specify:

 - The path to the directory containing the genome reference index files.
 - A file containing a set of SRA run IDs you want to download or the path were FASTQ files are stored on the local system.

 For example:

.. code:: bash

  nextflow run systemsgenetics/gemmaker -profile singularity \
    --pipeline salmon \
    --salmon_index_path Arabidopsis_thaliana.TAIR10.salmon.indexed \
    --sras SRAs.txt

Use Hisat2
..........
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

Use Local FASTQ Files
.....................
If your FASTQ files are local to your computer you must provide the ``--input`` argument when launching Nextflow and indicate the `GLOB pattern <https://en.wikipedia.org/wiki/Glob_(programming)>`_ than is needed to find the files:

.. code:: bash

  nextflow run systemsgenetics/gemmaker -profile singularity \
    --pipeline kallisto \
    --kallisto_index_path Arabidopsis_thaliana.TAIR10.kallisto.indexed \
    --input "../../01-input_data/RNA-seq/fastq/*{1,2}.fastq"

In the example above the ``--input`` argument indicates that FASTQ files are found in the ``../../01-input_data/RNA-seq/fastq/`` directory and GEMmaker should use all files that match the GLOB pattern ``*{1,2}.fastq``.

.. note ::

  GEMmaker currently expects that all fASTQ files have a `1` or `2` suffix. For paired files two files with the same name but each suffix respectively.

Use Both Local and SRA Files
............................
You can combine data from the NCBI SRA with local files in a single run of GEMmaker by providing both the ``--sras`` and ``--input`` arguments.

.. code:: bash

  nextflow run systemsgenetics/gemmaker -profile singularity \
    --pipeline kallisto \
    --kallisto_index_path Arabidopsis_thaliana.TAIR10.kallisto.indexed \
    --input "../../01-input_data/RNA-seq/fastq/*{1,2}.fastq" \
    --sras SRAs.txt

Using Paired-End Local Data
...........................
If your data is paired-end you must provide a `GLOB <https://en.wikipedia.org/wiki/Glob_(programming)>`_ pattern for the ``--input`` argument that can distinguish between the sample name and the suffix that indicates the pair.  Usually, paired-files have a ``1.fastq`` or ``2.fastq`` suffix on all file names.  Therefore, the GLOB given example given above is appropriate: ``*{1,2}.fastq``. The ``{1,2}`` indicates where the ``1`` and ``2`` are at in file name. However, if your files are named differently, be sure to use a GLOB pattern that can differentiate the pairs.

.. warning ::

    If the GLOB you provide cannot distinguish between pairs then GEMmaker will treat them as non-paired.

Using Non Paired-End Local Data
...............................
If your data is not paired-end then the `GLOB <https://en.wikipedia.org/wiki/Glob_(programming)>`_ pattern for the ``--input`` argument simply needs to find all of the FASTQ files.  For example, if your FASTQ files have a ``.fastq`` suffix the following GLOB would be appropriate:  ``*.fastq"``.

Using Both Paired-End and Non Paired Local Data
...............................................
GEMmaker can work with both paired and non-paired data in the same data set. The only stipulation is that the non-paired data must follow the same naming convention as the paired data. See the section `Using Paired-End Local data`_. For example, if your paired files have a ``1.fastq`` and ``2.fastq`` extension, then the non-paired files should have a ``1.fastq`` suffix as well.

Resuming After Failure
''''''''''''''''''''''
If for some reason GEMmaker fails to fully complete and Nextflow reports some form of error. You can resume execution of the workflow, afer correcting any problems, by passing the ``-resume`` flag to GEMmaker. For example to resume a failed Kallisto run:

.. code:: bash

  nextflow run systemsgenetics/gemmaker -profile singularity \
    -resume \
    --pipeline kallisto \
    --kallisto_index_path Arabidopsis_thaliana.TAIR10.kallisto.indexed \
    --sras SRAs.txt

GEMmaker should resume processing of samples without starting over.

Skipping Samples
................
You may find that a sample is problematic. It may be corrupt, does not align or has other problems that may cause GEMaker to fail. For such samples that cause GEMmaker to fail, you have two options. You can either remove the bad samples and restart GEMmaker or you can resume, as just described in the previous section, but first add the sample names to a new file, one per line, then, use the ``--skip_samples`` argument to tell GEMmaker about this file.  For example:

.. code:: bash

  nextflow run systemsgenetics/gemmaker -profile singularity \
    --pipeline kallisto \
    --kallisto_index_path Arabidopsis_thaliana.TAIR10.kallisto.indexed \
    --sras SRAs.txt \
    --skip_samples samples2skip.txt

In the example above any samples that should be skipped should be added to the ``samples2skip.txt`` file.

.. warning ::

    Note, when you provide SRA IDs to GEMmaker you provide the RUN IDs, but multiple run IDs can be contained in a single sample.  To skip a sample, you must provide the sample ID. For SRA, these  begin with the prefix SRX, DRX or ERX, where as run IDs begin with SRR, DRR or ERR.

Running on a Cluster
''''''''''''''''''''
If you want to run GEMmaker on a local High Performance Computing Cluster (HPC) that uses a scheduler such as SLURM or PBS, you must first create a configuration file to help GEMmaker know how to submit jobs.  The file should be named ``nextflow.config`` and be placed in the same directory where you are running GEMmaker.  Below is an example ``nextflow.config`` file for executing GEMmaker on a cluster that uses the SLURM scheduler.

.. code::

   profiles {
      my_cluster {
         process {
            executor = "slurm"
            queue = "<queue name>"
            clusterOptions = ""
         }
         executor {
            queueSize = 120
        }
      }
   }

In the example above we created a new profile named ``my_cluster``. Within the stanza, the placeholder text ``<queue name>`` should be replaced with the name of the queue on which you are allowed to submit jobs. If you need to provide specific options that you would normally provide in a SLURM submission script (such as an account or other node targetting settings) you can use the ``clusterOptions`` setting.

Next, is an example SLURM submission script for submitting a job to run GEMmaker. Please note, this is just an example and your specific cluster may require slightly different configuration/usage. The script assumes your cluster uses the lmod system for specifying software.

.. code:: bash

    #!/bin/sh
    #SBATCH --partition=<queue_name>
    #SBATCH --nodes=1
    #SBATCH --ntasks-per-node=1
    #SBATCH --time=10:00:00
    #SBATCH --job-name=GEMmaker
    #SBATCH --output=%x-%j.out
    #SBATCH --error=%x-%j.err

    module add java nextflow singularity

    nextflow run systemsgenetics/gemmaker \
      -profile my_cluster,singularity \
      -resume \
      --pipeline kallisto \
      --kallisto_index_path Araport11_genes.201606.cdna.indexed \
      --sras  SRA_IDs.txt \
      --max_cpus 120

Notice in the call to nextflow, the profile ``my_cluster`` has been added along with ``singularity``, also, the ``--max_cpus`` argument has been set to the same size as the ``queueSize`` value in the config file. The default value of ``--max_cpus`` is 4 and won't allow the workflow to expand beyond 4 CPUs if it is not increased to match the config file.


Intermediate Files
''''''''''''''''''
GEMmaker was designed to limit the storage requirements in order to allow for processing of large numbers of FASTQ files without overrunning storage requirement.  By default it will remove all large intermediate files to keep space usage to a minimum. However, you can indicate what intermediate files you would like to keep by providing any of the following arguments and setting them to ``true``.  For example, to keep the downloaded SRA files the ``keep_sra`` argument would be provided and set to true:

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


Configuration
'''''''''''''
The instructions above provide details for running GEMmaker using Singularity. For most instances you probably won't need to make customizations to the workflow configuration. However, should you need to, GEMmaker is a `nf-core <https://nf-co.re/>`_ compatible workflow.  Therefore, it follows the general approach for workflow configuration which is described at the `nf-core Pipeline Configuration page <https://nf-co.re/usage/configuration>`_.  Please see those instructions for the various platforms and settings you can configure.  However, below are some quick tips for tweaking GEMmaker.

In all cases, if you need to set some customizations you must first create a configuration file.  The file should be named ``nextflow.config`` and be placed in the same directory where you are running GEMmaker.

Configuration for a Cluster
...........................
To run GEMmaker on a computational cluster you will need to to create a custom configuration.  Instructions and examples are provided in the `Running on a Cluster`_ section.

Increasing Resources
.....................
You may find that default resources are not adequate for the size of your data set.  You can alter resources requested for each step of the GEMmaker workflow by using the ``withLabel`` scope selector in a custom ``nextflow.config`` file.

For example, if you have thousands of SRA data sets to process, you may need more memory allocated to the ``retrieve_sra_metadata`` step of the workflow. All steps in the workflow have a "label" that you can use to indicate which step resources should be changed. Below is an example ``nextflow.config`` file where a new profile named ``custom`` is provided where the memory has been increased for the ``retrieve_sra_metadata``.

.. code::

    profiles {
        custom {
            process {
                withLabel:retrieve_sra_metadata {
                    memory = "10.GB"
         	    }
            }
        }
    }

This new ``custom`` profile can be used when calling GEMmaker. The following is an example Kallisto run of GEMmaker using the custom and singularity profiles:

.. code:: bash

  nextflow run systemsgenetics/gemmaker -profile custom,singularity \
    --pipeline kallisto \
    --kallisto_index_path Arabidopsis_thaliana.TAIR10.kallisto.indexed \
    --sras SRAs.txt

Nextflow provides many "directives", such as ``memory`` that you can use to alter or customize the resources of any step (or process) in the workflow.  You can find more about these in the `Nextflow documentation. <https://www.nextflow.io/docs/latest/process.html#directives>`_ Some useful directives are:

- `memory <https://www.nextflow.io/docs/latest/process.html#memory>`_: change the amount of memory allocated to the step.
- `time <https://www.nextflow.io/docs/latest/process.html#time>`_: change the amount of time allocated to the step.
- `disk <https://www.nextflow.io/docs/latest/process.html#disk>`_: defines how much local storage is required.
- `cpus <https://www.nextflow.io/docs/latest/process.html#cpus>`_: defines how many threads (or CPUs) the task can use.

The "labels" that GEMmaker provides and which you can set custom directives include:

- ``retrieve_sra_metadata``:  For the step that retrieves metadata from the NCBI web services for the SRR run IDs that were provided. This step can require more memory than the defaults if there are huge numbers of samples.
- ``download_runs``: For the step is used for downloading SRA files from NCBI.
- ``fastq_dump``: For the step that is used after downloading SRA files and converting them to FASTQ files.
- ``fastqc``: For the step where the FastQC program is used which generates quality reports on FASTQ files.
- ``kallisto``: For the step the runs the Kallisto tool.
- ``salmon``: For the step that runs the Salmon tool.
- ``trimmomatic``: For the step that runs the Trimmomatic step which only runs when hisat2 is the desired pipeline.
- ``hisat2``: For the step that runs the hisat2 tool.
- ``samtools``: For the step that runs when the samtools tool is used after Hisat2 runs. This step only runs when the hisat2 pipeline is used.
- ``stringtie``: For the step that runs the stringtie tool and which only runs when the hisat2 pipeline is used.
- ``multiqc``: For the step that runs the MultiQC results summary report.
- ``create_gem``: For the step that creates the final GEM files.
- ``multithreaded``:  For all of the tools that support multithreading you can use this label to set a default number of CPUs using the ``cpus`` directive.  These tools include Salmon, Kallisto, Trimmomatic, Hisat2 and Stringtie.  By using this label you set set the same number of ``cpus`` for all multithreaded steps at once.

Using the Development Version
'''''''''''''''''''''''''''''
New updates to GEMmaker, prior to issuing a formal release, are held in the ``dev`` branch of the GEMmaker github repository. It is recommended to always use a formal release of GEMmaker, however, you can test the most recent improvements prior to release.  To do so, use the ``-r dev`` argument when running GEMmaker. For example:

.. code:: bash

  nextflow run systemsgenetics/gemmaker -r dev -profile singularity \
    --pipeline kallisto \
    --kallisto_index_path Arabidopsis_thaliana.TAIR10.kallisto.indexed \
    --sras SRAs.txt

The ``-r dev`` argument forces Nextflow to use the development version of GEMmaker rather than the most recent stable version.

.. note::

    You can find the most recent documentation for the ``dev`` branch at https://gemmaker.readthedocs.io/en/dev/
