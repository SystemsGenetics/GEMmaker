[![DOI](https://zenodo.org/badge/114067776.svg)](https://zenodo.org/badge/latestdoi/114067776)


![GEMmaker Logo](images/GEMmaker-logo-sm.png)


GEMmaker is a [Nextflow](https://www.nextflow.io/) workflow for large-scale gene expression sample processing, expression-level quantification and Gene Expression Matrix (GEM) construction. Results from GEMmaker are useful for differential gene expression (DGE) and gene co-expression network (GCN) analyses. This report is the MultiQC summary of the GEMmaker workflow results. The GEMmaker workflow currently supports Illumina RNA-seq datasets. 

The following flowchart describes the workflow that GEMmaker provides:


![flowchart](images/flowchartgen.png)


## Prerequisites

Before execution of GEMmaker you must have the necessary software. The following list provides the set of tools and versions that have been verified to work with GEMmaker. Note: newer versions of these tools are assumed to also work. Older versions may work but have not been tested:

- [Python3](https://www.python.org)
    - [xmltodict](https://github.com/martinblech/xmltodict) Python package used here for parsing of NCBI metadata
- [NextFlow](https://www.nextflow.io/) v0.28:  executes the workflow.
- [sratoolkit](https://www.ncbi.nlm.nih.gov/books/NBK158900/) v2.8.0:  Downloads SRA files from NCBI using the SRA Run IDs.
- [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) v0.11.7:  Generates read quality statistics for FASTQ files used by the workflow.
- [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic) v0.38:  Removes low-quality bases from the ends of reads and removes adapter sequences.
- [Hisat2](https://ccb.jhu.edu/software/hisat2/index.shtml) v2.1.0:  Aligns cleaned reads to the reference genome.
- [Samtools](http://www.htslib.org/) v1.3.1:  Used for indexing and sorting of BAM files created by Hisat2.
- [StringTie](http://www.ccb.jhu.edu/software/stringtie/) v1.3.4d:  Performs gene expression quantification.
- [MultiQC](http://multiqc.info/) (optional) v1.5:  If MultiQC is installed it will generate a full summary report for the entire workflow.

Additionally, GEMMaker requires that these software be available on your computer (or cluster) via the [Lmod](https://www.tacc.utexas.edu/research-development/tacc-projects/lmod) system and expects that the modules for these software are named as follows:

- sratoolkit
- fastQC
- trimmomatic
- hisat2
- stringtie
- sammtools

Preparing the prerequisites may be challenging for some.  Additionally, some HPC systems may have the software available but not using the module names listed above.  There are a few options to simplify execution of GEMmaker despite these problems

### Stand-alone computer

You can execute GEMmaker using a stand alone computer by using the [GEMmaker-docker image](https://github.com/SystemsGenetics/GEMmaker-docker).  See the GitHub site for instructions for execution of your workflow using that image.  It contains all of the necessary software needed to execute the workflow. No installation of software dependencies is reqiured and it will not conflict with existing software.  Note that execution of GEMmaker for a large number of samples on a single stand-alone machine is not recommended as it may take a very long time to complete.

### High-Performance Computing (HPC) cluster

To execute GEMmaker on an HPC cluster you must do **only one** of the following:
1) Ask your HPC admins to install the necessary software using the module names specified above.
2) Install the software into your own space and create your own module files that provide names for the software.  You can find examples of module files in the 'files' directory of the [GEMmaker-docker](https://github.com/SystemsGenetics/GEMmaker-docker) repository.
3) Edit GEMmaker's main.nf script and alter the module names to match those of your HPC system.

## Preparing the Workflow

After ensuring that all necessary software prerequisites are available, clone GEMmaker into a working directory.  

To clone the workflow into a directory:
```bash
nextflow clone SystemsGenetics/GEMmaker target-dir
```

As with all NextFlow workflows, you can configure the behavior of the workflow by creating a **nextflow.config** file.  The GEMmaker workflow provides an example file (`nextflow.config.example`) you can copy to get started. 
```bash
cp nextflow.config.example nextflow.config
```

Now, edit `nextflow.config` according to the inline instructions.  You may want to refer to the [Nextflow configuration documentation](https://www.nextflow.io/docs/latest/config.html) to set proper "profile" settings for your computing infrastructure. For example, to execute this workflow on an HPC system you must provide the name of the job scheduler and the job queue to submit to.

## Preparing Your Data

### Sample Files

GEMmaker supports processing of sample files that are already present on your local computer or samples that are already stored in the [NCBI SRA repository](https://www.ncbi.nlm.nih.gov/sra). 

- For local sample files, identify a [glob pattern](https://en.wikipedia.org/wiki/Glob_(programming)) that finds these files.  
- For samples on NCBI, identify the NCBI SRA run IDs of the samples you want to analyze. The Run IDs typically start with an SRR, ERR or DRR prefix. hese sample run IDs must be placed, one per line, in a file and the filename should be set in the **remote_list_path** of `nextflow.config`.

__NOTE__: The SRA Toolkit caches SRA files in your home directory by default. For large experiments this cache can become quite large, which may become an issue on some HPC clusters where each user is given a disk quota for their home directory. You can change the location of the cache by running `vdb-config -i` (see [SRA Toolkit Configuration](https://github.com/ncbi/sra-tools/wiki/Toolkit-Configuration)).

### Reference Genome Files

Download the genome annotation/reference files.  You must have the following:

- A FASTA file containing the full genomic sequence (either pseudomolecules or scaffolds). Note, if your genome file is extremely large with hundreds of thousands of contigs/scaffolds, you may want to reduce the size of the FASTA file to contain only those contigs/scaffolds with predicted annotated genes.
- A [GTF](https://uswest.ensembl.org/info/website/upload/gff.html) file containing the gene models. Sometimes a genome assembly does not provide a GTF file, but rather provides a [GFF3](https://uswest.ensembl.org/info/website/upload/gff.html) file. You can convert the GFF file to a GTF file using the **gffread** program of [cufflinks](http://cole-trapnell-lab.github.io/cufflinks/file_formats/), which you may have to download separately.  An example command-line to convert a GFF3 to GTF is ```gffread [gff_file] -T -o [gtf_file]``` where `[gff_file]` and `[gtf_file]` should be substituted for the names of your GFF3 and desired GTF file respectively.

Additional considerations:

- You must have hisat2 index files of your genome sequence. These are constructed by using the **hisat2-build** command.
- The GTF file and the hisat2 index files must have the same prefix and this prefix must be identified in `nextflow.config` using the **prefix** parameter for **hisat2-build**.
- All of the genome annotation files must be in a directory and this directory must be identified in `nextflow.config` using the **ref** > **path** paramter.

As an example for a proper setup, you will notice that the GEM-maker project contains an **examples** directory and within the **examples/reference** directory all of the files have the same prefix of **CORG**. This prefix is set for the **prefix** parameter in `nextflow.config.example`. The example directory also contains an REMOTE_IDS.txt file containing a list of SRA fastq_run_IDs.

Once your files are prepared, you can execute the workflow.

## Executing the Workflow

To execute the workflow on a local machine use this command:
```bash
nextflow run main.nf -profile standard
```

To resume a workflow in the event of a failure:
```bash
nextflow run main.nf -profile standard -resume
```

To execute the workflow and generate trace, timeline and execution reports use this command:
```bash
nextflow run main.nf -profile standard -with-report -with-timeline -with-trace
```

To execute the workflow on a high performance compute cluter you must edit `nextflow.config` and add an appropriate profile for your system. Please see the [Nextflow documentation](https://www.nextflow.io/docs/latest/config.html#config-profiles).  Then repeat any of the commands above changing the -profile argument to use the new profile.

## Generating a Summary Report

The [MultiQC](http://multiqc.info) tool can be used with GEMmaker to generate a summary report of results from Trimmomatic, Hisat2 and samtools.  This report allows you to explore the quality of the data, trimming and alignments.  To generate the report you must have [MultiQC installed](http://multiqc.info/docs/#installing-multiqc).  Once installed, you can generate the report with the following command inside of the GEMmaker directory where your workflow was executed:

```bash
multiqc .
```

## Generating the Gene Expression Matrix (GEM)

After GEMmaker completes, the results for all steps for each sample are stored in directories specific for each sample.  You can find a Gene Expression Vector (GEV) for each sample in the sample directory. the GEV will be the file with the extension .fpkm (for FPKM values) or .tpm (for TPM values). To compile all GEVs into a Gene Expression Matrix (GEM) you can use the `create_GEM.py` script found in the scripts directory.  To see help documentation for this script run the following:

```bash
python ./scripts/create_GEM.py -h
```

To create a GEM file after GEMmaker completes using the TPM values use the following:

```bash
python ./scripts/create_GEM.py --source ./ --type TPM --prefix my_project
```

Be sure to change `my_project` to a meaningful prefix for your project. Once completed the GEM file named `my_project.GEM.TPM.txt` will be created.

If you want to combine the results of multiple GEMmaker runs into a single GEM you can do so.  This may be useful if you split a set of input files into different GEMmaker runs and now you want to combine then.  To do so, you can provide to the `--source` argument the path to each GEMmaker directory.

Once completed, you will have a new file named GEM.txt inside the working directory.
