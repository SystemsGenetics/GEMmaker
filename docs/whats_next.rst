What to do with the GEM?
------------------------
The Gene Expression Matrix (GEM) created by GEMmaker can be used for either Differential Gene Expression (DGE) analyssi or Gene Co-expression Network (GCN) analysis.

DGE Analysis
''''''''''''

The raw GEM can be used for differential gene expression (DGE) analysis in `edgeR <https://bioconductor.org/packages/release/bioc/html/edgeR.html>`__ and `DESeq2 <https://bioconductor.org/packages/release/bioc/html/DESeq2.html>`__.

Network Analysis
''''''''''''''''

The GEM can be used to construct a gene co-expression network (GCN). The most common tool for construction of GCNs is `WGCNA <https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/>`__.  However, the developers of GEMmaker have also developed a new tool for constructing condition-specific GCNs called `KINC <https://github.com/SystemsGenetics/KINC>`__. It is a high-performance application that can construct networks using Pearson or Spearman for pairwise correlation, as well as Gaussian mixture models (GMMs) for pairwise clustering. KINC is a Qt/`ACE <https://github.com/SystemsGenetics/ACE>`__ application that is capable of running on CPUs and GPUs, which means that it can scale to larger workloads.

.. note::

  Prior to network construction it is recommended to normalized (such as with quantile normalization) and log-transform the GEM.  GEMmaker does not provide this functionality.
