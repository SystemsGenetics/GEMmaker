.. _troubleshooting:

Troubleshooting
---------------
When running GEMmaker you may encounter the following issues.

Cannot get property 'remote_list_path' on null object
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
**Problem**: I received the following error:

.. code:: bash

  ERROR ~ Cannot get property 'remote_list_path' on null object

**Solution**: You forgot to create ``nextflow.config`` from the example config:

.. code:: bash

  cp nextflow.config.example nextflow.config

255 ERROR
~~~~~~~~~
On local machines, you may encounter the following Singularity error:

.. code:: bash

  ERROR  : Failed to set loop flags on loop device: Resource temporarily unavailable
  ABORT  : Retval = 255

This is caused by Singularity attempting to access the same image with more than one thread. The first process to access the image will lock it until it is read into memory. This can be safely ignored, as GEMmaker will automatically retry the process.

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

Get Help or Suggest Improvements
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

If you have questions, comments, suggestions for improvement or require help with setup and execution of GEMmaker please consider posting to the `GEMmaker issue board <https://github.com/SystemsGenetics/GEMmaker/issues>`_ on Github.
