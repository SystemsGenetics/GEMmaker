This README.md conatins information on what is going on behind the scene with the example data. To run the example data, please reference the README.md on the front page of this github repository.

This README.md is useful for understanding what you will need for a run of your own by providing examples.

# Example Data

GEMmaker comes with two examples **Local Run Example** and **Remote Run Example**. These are representative of the two types of runs GEMmaker is capable of:
* **Local Run** - used when you already have the raw fastq files on your local machine or local computing cluster. If the data is not published, use this type of run.
* **Remote Run** - used when you wish to have GEMmaker download files from NCBI automatically.

---

# Local Run Example
A **local run** is when you already have the raw fastq files on your local machine or computing cluster.

The Local Run Example uses the imaginary organism "Cool Organism" (CORG) and a data set of 3 artificially made RNA-seq runs for CORG. CORG has a very small "genome" of only 1,953 nucleotides, 2 "chromosomes" and 3 "genes". The 3 genes are named "gene\_Alpha", "gene\_Beta" and "gene\_Zeta".

This data set simulates a local run because the RNA-seq data is already on your local machine or cluster. This data set is very small, so it can be run on a desktop machine in a short amount of time.

### Reference directory
The reference directory for the **Local Run Example** is located at:
```bash
GEMmaker/examples/LocalRunExample/reference
```
This directory contains the made up reference genome file (**CORG.fna**), gtf file (**CORG.gtf**), and hisat files (**CORG.?/ht2**). It also has a **COMMANDS.sh** file with the command used to generate the hisat files from the genome file. You can modify this command to make hisat files for your own project.

### Data directory
The example data directory for the **Local Run Example** is located at:  
```bash
/GEMmaker/examples/LocalRunExample/Data/
```
This directory contains 3 "RNA-seq" data sets which are examples of unpaired data, and are each in a directory of their own. The file format for these reads is "?\_1.fastq" where the "?" is replaced by the number of the sample. GEM-maker finds these files through the glob pattern assigned to the "local\_samples\_path" in the **nextflow.config** file.

### Running Local Run Example Data
If you have not done so already, you can execute the Local Run Example workflow using the information contained in the section **Test using the example data** in the Readme.md on the front page of this github repository.

### Files Output by Local Run Example
Once executed, the Local Run Example should output 3 directories. These will be titled Sample_1, Sample_2, and Sample_3. In a real run, GEMmaker will automatically combine files that have the same experiment number( \[SED\]RX0000000 ) but different run numbers ( \[SED\]RR0000000 ), so it is possible that the \[SED\]RX number contains multiple \[SED\]RR runs. However, in the the local example, this is not the case.

In each output directory you will find the following files:
- **fastq**   The fastq reads file for the experiment
- **fastqc**  6 or 12 files (depending on paired or unpaired data) from fastqc. Fastqc is set up to check files before and after trimmomatic
- **sam**  alignment file
- **bam**   binary alignment file
- **ga**  expression level transcript abundance
- **fpkm**  2 column version of **ga** file with only gene and fpkm value

#### After the workflow

The output of GEM-maker can be used for several different analysis. The FPKM files can be combined into an expression matrix and then visualized using a heatmap. The following heatmap is the Local Example's fpkm values divided by 1000 in heatmap form. We can see that gene_Zeta remained constant across all three samples, gene_Beta decreased, and gene_Alpha increased.

![heatmap](../images/heatmap.png)

---
# Remote Run Example
A **Remote Run** is used when you wish to have GEMmaker download files from NCBI automatically

The Remote Run Example uses data from an unidentified bacteria that is located on NCBI. The data set is unusually small, which makes it useful as an example.

### Reference directory
The reference directory for the Local Run Example is located at:
```bash
GEMmaker/examples/RemoteRunExample/reference
```
It contains the reference files just like the directory for the Local Run Example, but these files are for the unidentified bacteria rather than CORG.

### Input File
Unlike a Local Run, a Remote Run requires a list of [NCBI Sequence Read Archive (SRA) IDs](https://www.ncbi.nlm.nih.gov/sra). These are unique identifiers of a nucleotide projects that are housed on NCBI. The file should contain one SRA per line, and not contain blank space at the bottom. Here is an example file:
```bash
GEMmaker/examples/RempRunExample/SRA_IDs.txt
```
The above file will allow GEMmaker to find and download the read for the unknown bacteria with the SRA ID [SRR649944](https://www.ncbi.nlm.nih.gov/sra/SRR649944/)
### Running Local Run Example Data
To run the Remote Run Example, you must change some parameters in the **nextflow.config** file.

First, you must change the **params.software_params.hisat2.path** so that it points to the proper reference files directory. It should be set to:
```bash
path = "${PWD}/examples/RemoteExampleRun/reference/"
```
The parameter right below this is **params.software_params.hisat2.prefix**. This should be changed to match the prefix of the reference files. In this example, the prefix should be set to:
```bash
prefix = "GCA_002793175.1_ASM279317v1_genomic"
```

**params.input.remote_list_path** should be set to point to SRA_IDs.txt file:
```bash
remote_list_path = "${PWD}/examples/RemoteExampleRun/SRA_IDs.txt"
```  
**params.input.local_samples_path** should be set to indicate that we have no local samples:
```bash
local_samples_path = "none"
```
### Files Output by Local Run Example
This will output 1 directory titled "SRR649944", which is the name of the run GEMmaker was set to download. It will contain the standard set of output files.

### Running both Local and Remote at the same time
GEMmaker is capable of running both Local and Remote at the same time. Just indicate where your local samples and list of remote samples is located in the **nextflow.config**
