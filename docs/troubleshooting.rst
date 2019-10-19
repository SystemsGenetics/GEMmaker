.. _troubleshooting:

Troubleshooting
---------------

ERROR  : Unknown image format/type
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
If you encounter the following error:

.. code:: bash

  ERROR  : Unknown image format/type: <directory to a singularity image file>
  ABORT  : Retval = 255

Most likely, Singularity encountered some problem when retrieving and building the software images that GEMmaker uses.  The solution is to just resume the GEMmaker workflow and the problem will most likely resolve itself.

ERROR  : No valid /bin/sh in container
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
If you encounter the following error message:

.. code:: bash

 ERROR  : No valid /bin/sh in container
 ABORT  : Retval = 255

Then this can be caused by any of the following problems:

-  Nextflow requires Singularity version 2. Using Singularity verion 3, for instance, will cause this problem.  
-  Nextflow will automatically download and build the singularity images for you.  If your `umask` is not set to create files that are readable and executable then you can get this error.  Setting a umask such as, ``umask u=rwx,g=rx,o=rx``, prior to running Nextflow will ensure the images are readable and executable.
-  The ``singularity.cacheDir`` setting in the ``nextflow.config`` file indicates where Nextflow will store the downloaded singularity images.  If you have a typo, an incorrect path, or say a space at the end of the path, then this will result in the image not being found and cause this error.

ERROR  : Failed to set loop flags on loop device: Resource temporarily unavailable
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
On local machines, you may encounter the following Singularity error:

.. code:: bash

  ERROR  : Failed to set loop flags on loop device: Resource temporarily unavailable
  ABORT  : Retval = 255

This is caused by Singularity attempting to access the same image with more than one thread. The first process to access the image will lock it until it is read into memory. This can be safely ignored, as GEMmaker will automatically retry the process.

Exception in thread "main" java.lang.OutOfMemoryError: Java heap space
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
This error can occur when GEMmaker runs the FastQC process and there is insufficient RAM available.  If running on a local machine you can adjust the ``params.execution.queue_size`` in the configuration file to decrease the number of concurrent jobs that run at a given time.  If running on an HPC system, you can increase the amount of memory requested for the job by altering the setting in the ``withlabel:fastqc`` section of the configuration file.  For example, to set the memory to 2 GB: 

.. code:: bash

   withLabel:fastqc {
     container = "gemmaker/fastqc:0.11.7-1.1"
     time = "24h"
     memory = "2 GB"
   }

Why is it taking so long to pull a docker/singularity image?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
This is dependent directly on your internet speed. The first time GEMmaker is run, it must download all of the programs it needs to run. This means it may take a little while longer to run the first time it is run on your machine.

What are all these WARN: messages?!
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
When running GEMmaker, you may see messages like this:

.. code:: bash

  WARN: The channel `create` method is deprecated -- it will be removed in a future release
  WARN: The `close` operator is deprecated -- it will be removed in a future release

Nextflow, the language GEMmaker is based on, is undergoing some upgrades. GEMmaker is based on a stable version of Nextflow that does not use the currently in progress changes. These errors can be safely ignored, as they are just warnings that Nextflow is undergoing upgrades.

GEMmaker seems hung and does not complete
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
If GEMmaker seems to have processed all the samples provided, but does not move on to the ``create_gem`` step it may be hung.  Sometimes this can occur if you have tried to run GEMmaker multiple times but changed the list of samples between runs.  The best solution is to only run GEMmaker with one set of samples, and to create a new installation of GEMmaker for other samples.  However, if you do not want to lose results, you can try to run the following to clear out the GEMmaker batch directories:

.. code:: bash

  rm -rf work/GEMmaker/*
  
SLURM:  exceeded memory limit
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
If you are launching GEMmaker on an HPC system with the SLURM scheduler you can sometimes get an `exceeded memory limit` similar to the following: 

.. code:: bash

  slurmstepd: error: Job 12254566 exceeded memory limit (7871840 > 6553600), being killed
  
If you have a lot of samples, Nextflow may need more memory.  Increasing the amount of memory in your SLURM submission script will correct this problem.  Remember to restart GEMmaker with the ``-resume`` flag to have it continue where it left off.

  
Get Help or Suggest Improvements
--------------------------------

If you have questions, comments, suggestions for improvement or require help with setup and execution of GEMmaker please consider posting to the `GEMmaker issue board <https://github.com/SystemsGenetics/GEMmaker/issues>`_ on Github.
