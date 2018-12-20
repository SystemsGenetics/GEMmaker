[![DOI](https://zenodo.org/badge/114067776.svg)](https://zenodo.org/badge/latestdoi/114067776)


![GEMmaker Logo](images/GEMmaker-logo-sm.png)


GEMmaker is a [Nextflow](https://www.nextflow.io/) workflow for large-scale gene expression sample processing, expression-level quantification and Gene Expression Matrix (GEM) construction. Results from GEMmaker are useful for differential gene expression (DGE) and gene co-expression network (GCN) analyses. This report is the MultiQC summary of the GEMmaker workflow results. The GEMmaker workflow currently supports Illumina RNA-seq datasets.

The following flowchart describes the workflow that GEMmaker provides:


![flowchart](images/flowchartgen.png)


## Prerequisites

Before execution of GEMmaker you must have the necessary software. The following list provides the set of tools and versions that have been verified to work with GEMmaker. NOTE: newer versions of these tools are assumed to also work, and older versions may work but have not been tested:

- [python3](https://www.python.org)
- [nextflow](https://www.nextflow.io/) v0.32:  Executes the workflow.
- [sratoolkit](https://www.ncbi.nlm.nih.gov/books/NBK158900/) v2.8.0:  Downloads SRA files from NCBI using the SRA Run IDs.
- [fastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) v0.11.7:  Generates read quality statistics for FASTQ files used by the workflow.
- [trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic) v0.38:  Removes low-quality bases from the ends of reads and removes adapter sequences.
- [hisat2](https://ccb.jhu.edu/software/hisat2/index.shtml) v2.1.0:  Aligns cleaned reads to the reference genome.
- [samtools](http://www.htslib.org/) v1.3.1:  Used for indexing and sorting of BAM files created by Hisat2.
- [stringTie](http://www.ccb.jhu.edu/software/stringtie/) v1.3.4d:  Performs gene expression quantification.
- [MultiQC](http://multiqc.info/) (optional) v1.5:  Generate a full summary report for the entire workflow.

Additionally, GEMmaker requires that some of these packages be available on your system as [Environment Modules](http://modules.sourceforge.net/) or [Lmod](https://www.tacc.utexas.edu/research-development/tacc-projects/lmod) modules and that these modules are named as follows:

- fastQC
- hisat2
- sammtools
- sratoolkit
- stringtie
- trimmomatic

Preparing the prerequisites may be challenging for some. Additionally, some HPC systems may have the software available but not using the module names listed above. There are a few options to simplify execution of GEMmaker despite these problems.

### Local machine

You can execute GEMmaker out-of-the-box using the [GEMmaker-docker](https://github.com/SystemsGenetics/GEMmaker-docker) image. Refer to the Github repo for usage instructions. It contains all of the necessary software needed to execute the workflow. No installation of software dependencies is required and it will not conflict with existing software. Note that execution of GEMmaker for a large number of samples on a local machine is not recommended as it may take a very long time to complete.

### High-Performance Computing (HPC) cluster

To execute GEMmaker on an HPC cluster you must do __only one__ of the following:
- Have your HPC admins install the necessary software using the module names specified above.
- Install the software into your own space and create your own module files for the software. You can find examples of module files in the `files` directory of the [GEMmaker-docker](https://github.com/SystemsGenetics/GEMmaker-docker) repository.
- Edit `main.nf` and alter the module names to match those of your HPC system.

## Preparing the Workflow

After ensuring that all necessary software prerequisites are available, clone GEMmaker into a working directory.

To clone the workflow into a directory:
```bash
nextflow clone SystemsGenetics/GEMmaker target-dir
```

As with all Nextflow workflows, you can configure the behavior of the workflow by creating a `nextflow.config` file.  The GEMmaker workflow provides an example file (`nextflow.config.example`) which you can copy to get started.
```bash
cp nextflow.config.example nextflow.config
```

Edit `nextflow.config` according to the inline instructions. At the very least, you will need to modify the input parameters and the execution profile. You may want to refer to the [Nextflow configuration documentation](https://www.nextflow.io/docs/latest/config.html) to set proper profile settings for your environment. For example, to run the workflow on an HPC system you will have to specify the "executor" that corresponds to your system's scheduler (such as `pbs`, `slurm`, etc), as well as any other properties specific to your system, such as your job queue.

## Preparing Your Data

### Sample Files

GEMmaker supports processing of sample files that are already present on your local computer or samples that are already stored in the [NCBI SRA repository](https://www.ncbi.nlm.nih.gov/sra).

- For local sample files, identify a [glob pattern](https://en.wikipedia.org/wiki/Glob_(programming)) that finds these files.
- For samples on NCBI, identify the NCBI SRA run IDs of the samples you want to analyze. The Run IDs typically start with an SRR, ERR or DRR prefix. hese sample run IDs must be placed, one per line, in a file and the filename should be set in the `remote_list_path` of `nextflow.config`.

__NOTE__: The SRA Toolkit caches SRA files in your home directory by default. For large experiments this cache can become quite large, which may become an issue on some HPC clusters where each user is given a disk quota for their home directory. You can change the location of the cache by running `vdb-config -i` (see [SRA Toolkit Configuration](https://github.com/ncbi/sra-tools/wiki/Toolkit-Configuration)).

### Reference Genome Files

Download the genome annotation/reference files.  You must have the following:

- A FASTA file containing the full genomic sequence (either pseudomolecules or scaffolds). Note, if your genome file is extremely large with hundreds of thousands of contigs/scaffolds, you may want to reduce the size of the FASTA file to contain only those contigs/scaffolds with predicted annotated genes.
- A [GTF](https://uswest.ensembl.org/info/website/upload/gff.html) file containing the gene models. Sometimes a genome assembly does not provide a GTF file, but rather provides a [GFF3](https://uswest.ensembl.org/info/website/upload/gff.html) file. You can convert the GFF file to a GTF file using the `gffread` program of [cufflinks](http://cole-trapnell-lab.github.io/cufflinks/file_formats/), which you may have to download separately.  An example command-line to convert a GFF3 to GTF is `gffread [gff_file] -T -o [gtf_file]` where `[gff_file]` and `[gtf_file]` should be substituted for the names of your GFF3 and desired GTF file respectively.

Additional considerations:

- You must have hisat2 index files of your genome sequence. These are constructed by using the `hisat2-build` command.
- The GTF file and the hisat2 index files must have the same prefix and this prefix must be identified in `nextflow.config` using the `prefix` parameter for `hisat2-build`.
- All of the genome annotation files must be in a directory and this directory must be identified in `nextflow.config` using the `ref > path` paramter.

The GEMmaker repo contains an `examples` directory which contains several small example setups. The `LocalRunExample` contains a `reference` directory whose files have the same prefix of `CORG`. This prefix is set for the `prefix` parameter in `nextflow.config.example`. The `RemoteRunExample` also contains an `SRA_IDS.txt` file which contains a list of SRA fastq_run_IDs to download from NCBI.

## Executing the Workflow

To execute the workflow on a local machine:
```bash
nextflow run main.nf -profile standard
```

To resume a workflow in the event of a failure:
```bash
nextflow run main.nf -profile standard -resume
```

To execute the workflow and generate trace, timeline and execution reports:
```bash
nextflow run main.nf -profile standard -with-report -with-timeline -with-trace
```

To execute the workflow on an HPC system you must edit `nextflow.config` and add an appropriate profile for your system. Refer to the [Nextflow documentation](https://www.nextflow.io/docs/latest/config.html#config-profiles). You can then use any of the above commands by changing the `-profile` argument to use your profile.

### Performance Considerations

For large experiments on an HPC system, it is important to make sure that you are effectively utilizing the resources of the system. There are a number of parameters in `nextflow.config` which can be used to increase performance based on the capabilities of your system:
- `params.execution.threads`: All processes which support multithreading (such as trimmomatic) will use this number of threads. This setting should be determined by the number of cores per node on your system; for example, if your system has nodes with 16 cores per node then you could set the number of threads to 16 to make full use of those nodes.
- `params.execution.queue_size`: Nextflow will only run 100 processes at a time by default, but you may be able to increase this value based on the queue limits of your system.
- `scratch`: If this directive is enabled in your profile, each process will attempt to use its own local disk space instead of the main working directory. If you are running the workflow from NFS storage then you may benefit greatly from this option, if the nodes on your system have sufficient local storage. This option can significantly reduce the amount of disk I/O on your NFS storage, since each process will only interact with the NFS storage at the beginning and end to transfer input and output files.

## Generating a Summary Report

The [MultiQC](http://multiqc.info) tool can be used with GEMmaker to generate a summary report of results from Trimmomatic, Hisat2 and samtools.  This report allows you to explore the quality of the data, trimming and alignments.  To generate the report you must have [MultiQC installed](http://multiqc.info/docs/#installing-multiqc).  Once installed, you can generate the report with the following command inside of the GEMmaker directory where your workflow was executed:

```bash
multiqc .
```

## Generating the Gene Expression Matrix (GEM)

After GEMmaker completes, the results for each sample are stored in a directory specific to that sample. The final output for each sample is a Gene Expression Vector (GEV) in the form of an FPKM or TPM file. To compile all GEVs into a Gene Expression Matrix (GEM) you can use the `create_GEM.py` script in the `scripts` directory.

To see help documentation for this script:
```bash
python ./scripts/create_GEM.py -h
```

To create a GEM file from the TPM files produced by GEMmaker:
```bash
python ./scripts/create_GEM.py --source ./ --type TPM --prefix my_project
```

The script will produce a GEM file called `my_project.GEM.TPM.txt`. Be sure to change `my_project` to a meaningful prefix for your project.

You can combine the results of multiple GEMmaker runs into a single GEM by providing a list of directories to the `--source` argument. This feature may be useful if you split a set of input files into several GEMmaker runs and now you need to combine then. The script will produce a file named `GEM.txt` in the working directory.



## Running with the Dockerfile

```bash
nextflow run main.nf -profile inDocker
```
