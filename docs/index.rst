
.. figure:: images/GEMmaker-logo-sm.png
   :alt: GEMmaker Logo

.. image:: https://zenodo.org/badge/114067776.svg
   :target: https://zenodo.org/badge/latestdoi/114067776

Welcome to GEMmaker's documentation!
====================================

GEMmaker is a `Nextflow <https://www.nextflow.io/>`_ workflow for large-scale gene expression sample processing, expression-level quantification and Gene Expression Matrix (GEM) construction. Results from GEMmaker are useful for differential gene expression (DGE) and gene co-expression network (GCN) analyses. The GEMmaker workflow currently supports Illumina RNA-seq datasets.

nf-core Compatibility
---------------------
GEMmaker is an `nf-core <https://nf-co.re/>`_ compatible workflow, however, GEMmaker is not an official nf-core workflow.  This is because nf-core offers the `nf-core/rnaseq <https://nf-co.re/rnaseq>`_ workflow which is an excellent workflow for RNA-seq analysis that provides similar functionality to GEMmaker.  However, GEMmaker is different in that it can scale to thousands of samples without exceeding local storage resources by running samples in batches and removing intermediate files.  It can do the same for smaller sample sets on machines with less computational resources.  This ability to scale is a unique adaption that is currently not provided by Nextflow.   When Nextflow does provide support for batching and scaling, the `nf-core/rnaseq <https://nf-co.re/rnaseq>`_ will be updated and GEMmaker will probably be retired in favor of the nf-core workflow. Until then, if you are limited by storage GEMmaker can help!

Acknowledgments
---------------
Development of GEMmaker was funded by the U.S. National Science Foundation Award `#1659300 <https://www.nsf.gov/awardsearch/showAward?AWD_ID=1659300&HistoricalAwards=false>`_.

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   installation
   explore
   prepare_reference
   prepare_samples
   execution
   results
   whats_next
   troubleshooting
