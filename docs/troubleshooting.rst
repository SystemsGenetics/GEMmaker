.. _troubleshooting:

Troubleshooting
---------------

GEMmaker is a very powerful tool, especially when coupled with its underlying technologes such as Nextflow, Docker/Singularity and HPC environments. But there are also many moving parts, and as a result there are a few things that can go wrong when you set up your experiments. These issues usually have simple solutions, so some of the more common issues are documented here.

**Problem**: I received the following error:

.. code:: bash

    ERROR ~ Cannot get property 'remote_list_path' on null object

**Solution**: You forgot to create `nextflow.config` from the example config:

.. code:: bash

    cp nextflow.config.example nextflow.config

Get Help or Suggest Improvements
--------------------------------

If you have questions, comments, suggestions for improvement or require help with setup and execution of GEMmaker please consider posting to the `GEMmaker issue queue <https://github.com/SystemsGenetics/GEMmaker/issues>`_.
