.. _execution:

Workflow Execution
------------------

Specifying an Environment
`````````````````````````

GEMmaker is designed to run on a number of different environments, including:

- Local machine
- HPC cluster
- Kubernetes cluster

Additionally, each of these environments may be able to provide software dependencies through Docker, Singularity, or Environment Modules. Here we provide some example command lines for each of these scenarios.

**Local Machine**

To run GEMmaker on a local machine:

.. code:: bash

  # assume software dependencies are installed locally ...
  nextflow run main.nf -profile standard

  # ... or use docker
  nextflow run main.nf -profile standard,docker

**HPC Cluster**

To execute the workflow on an HPC system you will likely need to edit ``nextflow.config`` and add an appropriate profile for your system. Refer to the `Nextflow documentation <https://www.nextflow.io/docs/latest/executor.html>`__ for the list of available executors. Here we provide some examples based on the existing profiles.

To run GEMmaker on Kamiak:

.. code:: bash

  # use environment modules ...
  nextflow run main.nf -profile slurm,modules_kamiak

  # ... or singularity
  nextflow run main.nf -profile slurm,singularity

To run GEMmaker on Palmetto:

.. code:: bash

  # use environment modules ...
  nextflow run main.nf -profile pbs,modules_palmetto

  # ... or singularity
  nextflow run main.nf -profile pbs,singularity

**Kubernetes Cluster**

GEMmaker can be run on a `Kubernetes <https://kubernetes.io/>`__ cluster with minimal effort, but there are additional steps required to configure the cluster and transfer input data and output data before and after execution. Consult the `kube-runner <https://github.com/SystemsGenetics/kube-runner>`__ project for instructions.

Performance Considerations
==========================

For large experiments on an HPC system, it is important to make sure that you are effectively utilizing the resources of the system. There are a few settings in ``nextflow.config`` which can be used to maximize performance based on the capabilities of your system:

- **Multithreading**: Processes which support multithreading (such as trimmomatic) will use multiple threads according to the number of CPUs allocated to the process. Refer to the ``pbs`` and ``slurm`` profiles for examples of how to allocate more CPUs for multithreaded processes. This setting should be determined by the number of cores per node on your system; for example, if your system has nodes with 16 cores per node then you could set the number of threads to 16 to make full use of those nodes. Note, however, that you may also need to consider the memory available on each node, as well as the potentially higher queueing time for jobs that request more resources.

- **Queue size**: Nextflow will only run up to 100 processes at a time by default (``params.execution.queue_size``), but you may be able to increase this value based on the queue limits of your system.

Resuming a Previous Run
=======================

In the event of a failure you can resume a previous workflow run:

.. code:: bash

  nextflow run main.nf -resume

Generating a Summary Report
===========================

The `MultiQC <http://multiqc.info>`__ tool will automatically generate a report on how each process ran.

Generating a Gene Expression Matrix (GEM)
=========================================

After GEMmaker completes, the resulting GEMs will be output to ``output/GEMs/`` by default. This directory contains the final gene-expression matrices in raw, TPM and FPKM form, depending on which output formats are enabled in ``nextflow.config``.

Using GEMs in Other Workflows
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

DGE Analysis
============

The raw GEM can be used for DGE analysis in edgeR and other DGE software.

Network Analysis
================

Any GEM can be used to construct a gene-coexpression network (GCN). `KINC <https://github.com/SystemsGenetics/KINC>`__ is a high-performance application that can construct networks using Pearson or Spearman for pairwise correlation, as well as Gassian mixture models (GMMs) for pairwise clustering. KINC is a Qt/`ACE <https://github.com/SystemsGenetics/ACE>`__ application that is capable of running on CPUs and GPUs, which means that it can scale to larger workloads.
