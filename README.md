[![DOI](https://zenodo.org/badge/114067776.svg)](https://zenodo.org/badge/latestdoi/114067776)


![GEMmaker Logo](images/GEMmaker-logo-sm.png)


The GEMmaker project is a [NextFlow](https://www.nextflow.io/) workflow that generates a file containing FPKM values for all genes in each sample of an RNA-seq data set.
In other words, a Gene Expression Vector (GEV) is created for each sample. GEMmaker can automatically download samples from [NCBI SRA](https://www.ncbi.nlm.nih.gov/sra), or can operate on files that are stored locally.
This workflow combines the [sratoolkit](https://www.ncbi.nlm.nih.gov/books/NBK158900/), [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/), [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic), [Hisat2](https://ccb.jhu.edu/software/hisat2/index.shtml), [Samtools](http://www.htslib.org/), and [StringTie](http://www.ccb.jhu.edu/software/stringtie/) software packages.
The workflow expects the Lua-based [Lmod](https://lmod.readthedocs.io/en/latest/) software module system is installed with each software described above, making them available via the module system. In addition, user must have python3 and the python package [xmltodict](https://github.com/martinblech/xmltodict) installed.
The GEMmaker workflow currently supports Illumina RNA-seq data sets. It is intended to be run on a high-performance compute cluster.

For testing purpose, or for execution of a small data set (or large data set if sufficient storage is available), a Docker image is available that contains all of the necessary software components: https://github.com/SystemsGenetics/GEMmaker-docker

Note: The GEMmaker worflow is not configured to use Hisat2/Stringtie to identify novel splice varients or gene models.
It uses the existing predicted gene models as provided by a reference genome's assembly annotation.  The following flowchart describes the steps in this workflow:

![flowchart](images/flowchartgen.png)

---


## Prepare the Workflow

First, clone this workflow project into a working directory.  

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

GEMmaker comes with two examples **Local Run Example** and **Remote Run Example**

To execute the GEMmaker with an example data set you must first rename the **nextflow.config.example** file as **nextflow.config**.

You should then ensure that the **trimmomatic.clip_path** option in the **nextflow.config** file is set to the full path where the Trimmomatic clipping files are housed.  Replace the text **<ILLUMINACLIP_PATH>** placeholder text with the path.

The example config file also has an example profile for running this workflow on a SLURM cluster. If you want to use the SLURM profile you must, at a minimum, change the **<QUEUE_NAME>** placeholder text to be the name of the queue used for submission on your cluster.  If you require additional settings you can adjust the profile as per the [NextFlow configuration documentation](https://www.nextflow.io/docs/latest/config.html#config-profiles).

Additional information on what is in the test data can be found in the **GEMmaker/examples/README.md** file. This contains information on what is going on behind the scene, and what to do if you are interested in using the **Remote Run** feature of GEMmaker, which will automatically download RNA-seq read files from NCBI.
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
