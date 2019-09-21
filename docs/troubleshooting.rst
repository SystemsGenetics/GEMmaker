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


Get Help or Suggest Improvements
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

If you have questions, comments, suggestions for improvement or require help with setup and execution of GEMmaker please consider posting to the `GEMmaker issue board <https://github.com/SystemsGenetics/GEMmaker/issues>`_ on Github.
