[![DOI](https://zenodo.org/badge/114067776.svg)](https://zenodo.org/badge/latestdoi/114067776)


![GEMmaker Logo](images/GEMmaker-logo-sm.png)


GEMmaker is a [Nextflow](https://www.nextflow.io/) workflow for large-scale gene expression sample processing, expression-level quantification and Gene Expression Matrix (GEM) construction. Results from GEMmaker are useful for differential gene expression (DGE) and gene co-expression network (GCN) analyses. This report is the MultiQC summary of the GEMmaker workflow results. The GEMmaker workflow currently supports Illumina RNA-seq datasets. 

The following flowchart describes the workflow that GEMmaker provides:

![flowchart](images/flowchartgen.png)

---

## Prerequisites
Before execution of GEMmaker you must have the necessary software. The following list provides the set of tools and versions that have been verified to work with GEMmaker. Note: newer versions of these tools are assumed to also work. Older versions may work but have not been tested:

- [Python3](https://www.python.org)
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

#### Stand-alone computer
You can execute GEMmaker using a stand alone computer by using the [GEMmaker-docker image](https://github.com/SystemsGenetics/GEMmaker-docker).  See the GitHub site for instructions for execution of your workflow using that image.  It contains all of the necessary software needed to execute the workflow. No installation of software dependencies is reqiured and it will not conflict with existing software.  Note that execution of GEMmaker for a large number of samples on a single stand-alone machine is not recommended as it may take a very long time to complete.

#### High-Performance Computing (HPC) cluster
To execute GEMmaker on an HPC cluster you must do one of the following:
1) Ask your HPC admins to install the necessary software using the module names specified above.
2) Install the softwqare into your own space and create your own module files that provide names for the software.  You can find examples of module files in the 'files' directory of the [GEMmaker-docker](https://github.com/SystemsGenetics/GEMmaker-docker) repository.
3) Edit the main.nf script an alter the module names to match those of your HPC system.

## Prepare the Workflow

After ensuring that all necessary software prerequisites are available, clone GEMmaker into a working directory.  

To clone the workflow into a directory:
```bash
nextflow clone SystemsGenetics/GEMmaker target-dir
```
As with all NextFlow workflows, you can configure the behavior of the workflow by creating a **nextflow.config** file.  The GEMmaker workflow provides an example file (nextflow.config.example) you can rename to get started. 

```bash
mv nextflow.config.example nextflow.config
```
Now, edit the nextflow.config file according to the inline instructions and the [Nextflow configuration documentation](https://www.nextflow.io/docs/latest/config.html)

---

## Test using the example data

GEMmaker comes with two examples **Local Example** and **Remote Example**

To execute the GEMmaker with an example dataset you must first rename the **nextflow.config.example** file as **nextflow.config**.

You should then ensure that the **trimmomatic.clip_path** option in the **nextflow.config** file is set to the full path where the Trimmomatic clipping files are housed.  Replace the text **<ILLUMINACLIP_PATH>** placeholder text with the path.

The example config file also has an example profile for running this workflow on a SLURM cluster. To use the SLURM profile you must, at a minimum, change the **<QUEUE_NAME>** placeholder text to be the name of the queue used for submission on your cluster.  If you require additional settings you can adjust the profile as per the [NextFlow configuration documentation](https://www.nextflow.io/docs/latest/config.html#config-profiles).

### Local Example

GEM-maker comes with a Local example dataset, which is stored in the **examples** folder. The **nextflow.config.example** (which should be renamed **nextflow.config** to run the example) is set up to run this data when you clone the repository.

 The data is from an imaginary organism with the name of "Cool Organism", which is abbreviated "CORG". CORG has a very small "genome" of only 1471 nucleotides, 2 "chromosomes" and 3 "genes". The 3 genes are named "gene\_Alpha", "gene\_Beta" and "gene\_Zeta". The made up reference genome file **CORG.fna**, gtf file **CORG.gtf**, and hisat files **CORG.?/ht2** are stored in the directory "./GEM-maker/examples/reference/".

The example data consists of 3 "RNA-seq" data sets which are contained in the directory "./GEM-maker/examples/Data/". They are examples of unpaired data, and are each in a folder of their own. The file format for these reads is "?\_1.fastq" where the "?" is replaced by the number of the sample. GEM-maker finds these files through the glob pattern assigned to the "local\_samples\_path" in the **nextflow.config** file.

Once you understand the above information, run the Local example dataset using the commands in the section below that is titled **Executing the Workflow**

Once executed, the local example should output 3 directories. GEM-maker will automatically combine files that have the same experiment number( \[SED\]RX0000000 ) but different run numbers ( \[SED\]RR0000000 ), so it is possible that the \[SED\]RX number contains multiple \[SED\]RR runs. However, in the the local example, this is not the case.

In each output directory you will find the following files:
- **fastq**   The fastq reads file for the experiment
- **fastqc**  6 or 12 files (depending on paired or unpaired data) from fastqc. Fastqc is set up to check files before and after trimmomatic
- **sam**  alignment file
- **bam**   binary alignment file
- **ga**  expression level transcript abundance
- **fpkm**  2 column version of **ga** file with only gene and fpkm value

#### After the workflow

The output of GEM-maker can be used for several different analysis. The FPKM files can be combined into an expression matrix and then visualized using a heatmap. The following heatmap is the Local Example's fpkm values divided by 1000 in heatmap form. We can see that gene_Zeta remained constant across all three samples, gene_Beta decreased, and gene_Alpha increased.

![heatmap](images/heatmap.png)

### Remote Example
If you wish to use GEMmaker to download all or some of your fastq files from NCBI, you would also need to include an REMOTE\_IDs.txt file. An example of such a file is located in "./GEMmaker/examples/". You can point GEMmaker to this list by modifying the "remote\_list\_path" parameter in the "nextflow.config" file.

See the **README.md** in the "/examples/" directory for more information about running the remote example.

---

## Test using your own data.

To prepare your own samples for execution you must peform the following:

- As with the example data set described above, you must edit the **netxflow.config** file and set the **trimmomatic.clip_path** and customize it for execution on a cluster if desired.

- For local files, identify a glob pattern that you can add to the **nextflow.config** file that will identify these files.

- For files on NCBI, Identify in NCBI SRA the fastq_run_ids of the SRA samples you want to analyze.
  fastq_run_id numbers typically start with an SRR, ERR or DRR prefix.
  These sample run IDs must be placed, one per line, in a REMOTE_IDS.txt file.
  These will be downloaded automatically by the program.

- Download the genome annotation/reference files.
  You must have the following:

  - A FASTA file containing the full genomic sequence (either pseudomolecules or scaffolds). Note, if your genome file is extremely large with hundreds of thousands of contigs/scaffolds, you may want to reduce the size of the FASTA file to contain only those contigs/scaffolds with predicted annotated genes.

  - A GTF file containing the gene models. Sometimes a genome assembly does not provide a GTF file, but rather provides a GFF3 file. You can convert the GFF file to a GTF file using the **gffread** program of [cufflinks](http://cole-trapnell-lab.github.io/cufflinks/file_formats/), which you may have to download separately.  An example command-line to convert a GFF3 to GTF is ```gffread [gff_file] -T -o [gtf_file]``` where [gff_file] and [gtf_file] should be substituted for the names of your GFF3 and desired GTF file respectively.

  - You must have hisat2 index files of your genome sequence.
    These are constructed by using the **hast2-build** command.

  - Make sure that your GTF file has the exact same prefix as the hisat2 index files.

  - All of the genome annotation files must be in a directory and this directory must be identified in the **nextflow.config** file using the **ref** > **path** paramter.

- Finally, edit the **nextflow.config** file and change the **prefix** parameter to be the prefix used with **hisat2-build** when you created the index files.

As an example for a proper setup, you will notice that the GEM-maker project contains an **examples** directory and within the **examples/reference** directory all of the files have the same prefix of **CORG**.
You'll also notice that this same value is set for the **prefix** in the nextflow.config.example file.
The example directory also contains an REMOTE_IDS.txt file containing a list of SRA fastq_run_IDs.

Once your files are prepared, you can execute the workflow.

---

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
nextflow run main.nf -profile standard -resume -with-report execute-report.html -with-timeline timeline-report.html -with-trace
```

To execute the workflow on a high performance compute cluter you must edit the nextflow.config file and add an appropriate profile for your system. Please see the [Nextflow documentation] (https://www.nextflow.io/docs/latest/config.html#config-profiles).  Then repeat any of the commands above changing the -profile argument to use the new profile.

## Generating a Summary Report
The [MultiQC] (http://multiqc.info) tool can be used with GEMmaker to generate a summary report of results from Trimmomatic, Hisat2 and samtools.  This report allows you to explore the quality of the data, trimming and alignments.  To generate the report you must have [MultiQC installed] (http://multiqc.info/docs/#installing-multiqc).  Once installed, you can generate the report with the following command inside of the GEMmaker directory where your workflow was executed:

```bash
multiqc . 
```

## Generating the Gene Expression Matrix (GEM)
After GEMmaker completes, the results for all steps for each sample are stored in directories specific for each sample.  You can find a Gene Expression Vector (GEV) for each sample in the sample directory. the GEV will be the file with the .fpkm extension and contains the full vector of expression for all genes in genome.   To compile all of these GEVs into a Gene Expression Matrix execute the following script inside of the GEMmaker directory where your workflow was executed:

```bash
./scripts/fpkm2gem.sh
```
Once completed, you will have a new file named GEM.txt inside the working directory.
