.. _execution:

Step 5: Run GEMmaker
--------------------

Specifying an Environment
`````````````````````````

GEMmaker is designed to run on a number of different environments, including:

- Local machine
- HPC cluster
- Kubernetes cluster

Additionally, each of these environments may be able to provide software dependencies through Docker, Singularity, or Environment Modules. Here we provide some example command lines for each of these scenarios.

Local Machine
'''''''''''''
To run GEMmaker on a local machine:

.. code:: bash

  # assume software dependencies are installed locally ...
  nextflow run main.nf -profile standard

  # ... or use docker
  nextflow run main.nf -profile standard,docker

HPC Cluster
'''''''''''

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


Resuming a Previous Run
```````````````````````

In the event of a failure you can resume a previous workflow run:

.. code:: bash

  nextflow run main.nf -resume
