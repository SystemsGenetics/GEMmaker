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

prefetch issues
~~~~~~~~~~~~~~~
**Problem**: I received one of the following errors:

.. code:: bash

   transfer incomplete while reading file within network system module - Cannot KStreamRead:
   
or

.. code:: bash
  
   timeout exhausted while reading file within network system module - Cannot KStreamRead:

**Solution**
Check with your network engineers to explore potential points of slowness in your network infrastructure.  These errors probably occur when the network is too busy to pull files from NCBI SRA in a reasonable time.   Alternatively, edit the `nextflow.config` file and change the `maxRetries` setting to a higher level. If the problem is due to network slowness, attempting to retry may provide enough opportunities for a transfer to complete.

Get Help or Suggest Improvements
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

If you have questions, comments, suggestions for improvement or require help with setup and execution of GEMmaker please consider posting to the `GEMmaker issue board <https://github.com/SystemsGenetics/GEMmaker/issues>`_ on Github.
