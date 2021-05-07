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

-  Nextflow will automatically download and build the singularity images for you.  If your `umask` is not set to create files that are readable and executable then you can get this error.  Setting a umask such as, ``umask u=rwx,g=rx,o=rx``, prior to running Nextflow will ensure the images are readable and executable.


ERROR  : Failed to set loop flags on loop device: Resource temporarily unavailable
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
On local machines, you may encounter the following Singularity error:

.. code:: bash

  ERROR  : Failed to set loop flags on loop device: Resource temporarily unavailable
  ABORT  : Retval = 255

This is caused by Singularity attempting to access the same image with more than one thread. The first process to access the image will lock it until it is read into memory. This can be safely ignored, as GEMmaker will automatically retry the process.

Exception in thread "main" java.lang.OutOfMemoryError: Java heap space
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
This error can occur when GEMmaker runs the FastQC process and there is insufficient RAM available.  If running on a local machine you can adjust the ``--max_cpus`` argument to decrease the number of concurrent jobs that run at a given time.  If running on an HPC system, you can increase the amount of memory requested for the job by altering the setting in the ``--max_memory`` argument. See the `nf-core Pipeline Configuration page <https://nf-co.re/usage/configuration>`_ for more in-depth details for providing more custom configuration files.


Why is it taking so long to pull a docker/singularity image?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
This is dependent directly on your internet speed. The first time GEMmaker is run, it must download the Docker image it needs to run. This means it may take a little while longer to run the first time it is run on your machine.

GEMmaker seems hung and does not complete
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
If GEMmaker seems to have processed all the samples provided, but does not move on to the ``create_gem`` step it may be hung.  Sometimes this can occur if you have tried to run GEMmaker multiple times but changed the list of samples between runs.  The best solution is to only run GEMmaker with one set of samples in a single directory, and to use a different directory for other samples.

SLURM:  exceeded memory limit
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
If you are launching GEMmaker on an HPC system with the SLURM scheduler you can sometimes get an `exceeded memory limit` similar to the following:

.. code:: bash

  slurmstepd: error: Job 12254566 exceeded memory limit (7871840 > 6553600), being killed

If you have a lot of samples, Nextflow may need more memory.  Increasing the amount of memory in your SLURM submission script will correct this problem.  Remember to restart GEMmaker with the ``-resume`` flag to have it continue where it left off.


Get Help or Suggest Improvements
--------------------------------

If you have questions, comments, suggestions for improvement or require help with setup and execution of GEMmaker please consider posting to the `GEMmaker issue board <https://github.com/SystemsGenetics/GEMmaker/issues>`_ on Github.
