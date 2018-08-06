# Example Data

GEMmaker comes with two examples **Local Example Run** and **Remote Example Run**. Use a **Local Run** when you already have the raw fastq files on your local machine or local computing cluster. A **Remote Run** is when you wish to download files from NCBI automatically.

## Local Example Run
A **local run** is when you already have the raw fastq files on your local machine or computing cluster.

This directory contains both the reference material for the imaginary organism "Cool Organism" (CORG) and a dataset of 3 artificially made RNA-seq runs for CORG. This is simulating a local run becasue the RNA-seq data is already on your local machine or cluster. This data set is very small, so it can be run on a desktop machine with ease.

The reference directory contains the made up reference genome file (**CORG.fna**), gtf file (**CORG.gtf**), and hisat files (**CORG.?/ht2**). It also has a **COMMANDS.sh** file with the command used to generate the hisat files from the genome file.

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

![heatmap](../images/heatmap.png)

## Remote Example Run

If you wish to use GEMmaker to download all or some of your fastq files from NCBI, you would also need to include an REMOTE\_IDs.txt file. An example of such a file is located in "./GEMmaker/examples/". You can point GEMmaker to this list by modifying the "remote\_list\_path" parameter in the "nextflow.config" file.

This directory also contains a file titled "REMOTE_IDs.txt". This file is an example of a file that is in the proper format to download files from NCBI. To run this remote example, simpley change the "nextflow.config" file's "remote_list_path" parameter to be
```
'remote_list_path = "${PWD}/examples/REMOTE_IDs.txt"'
```
You can also change the "local_samples_path" to be "none" if you do not want the local files to run. This workflow can run both local and remote at the same time, however.

It should be noted that the output of the remote example is not as exciting as the local example, and is intended to be a demonstration of the download process works rather than the output process. The remote example uses the CORG reference files even though the remote data is *Arabidopsis thaliana*(Mouse Ear Crest). The final fpkm files will have values of 0 because of this, but the workflow will run as normal so that you can see how the remote process works. The reason that we have not included the *A. thaliana* reference files is to keep this git repository small.

It may be good practice for a new user of this workflow to generate the *A. thaliana* reference files that are associated with these small rna-seq datasets on their own. TAIR [(The Arabidopsis Information Resource)](https://www.arabidopsis.org/) and [ARAPORT](https://apps.araport.org/thalemine/begin.do) contain the files necessary to do this.
To execute the GEMmaker with an example dataset you must first rename the **nextflow.config.example** file as **nextflow.config**.

You should then ensure that the **trimmomatic.clip_path** option in the **nextflow.config** file is set to the full path where the Trimmomatic clipping files are housed.  Replace the text **<ILLUMINACLIP_PATH>** placeholder text with the path.

The example config file also has an example profile for running this workflow on a SLURM cluster. To use the SLURM profile you must, at a minimum, change the **<QUEUE_NAME>** placeholder text to be the name of the queue used for submission on your cluster.  If you require additional settings you can adjust the profile as per the [NextFlow configuration documentation](https://www.nextflow.io/docs/latest/config.html#config-profiles).
This directory also contains a file titled "REMOTE_IDs.txt". This file is an example of a file that is in the proper format to download files from NCBI. To run this remote example, simpley change the "nextflow.config" file's "remote_list_path" parameter to be
```
'remote_list_path = "${PWD}/examples/REMOTE_IDs.txt"'
```
You can also change the "local_samples_path" to be "none" if you do not want the local files to run. This workflow can run both local and remote at the same time, however.

It should be noted that the output of the remote example is not as exciting as the local example, and is intended to be a demonstration of the download process works rather than the output process. The remote example uses the CORG reference files even though the remote data is *Arabidopsis thaliana*(Mouse Ear Crest). The final fpkm files will have values of 0 because of this, but the workflow will run as normal so that you can see how the remote process works. The reason that we have not included the *A. thaliana* reference files is to keep this git repository small.

It may be good practice for a new user of this workflow to generate the *A. thaliana* reference files that are associated with these small rna-seq datasets on their own. TAIR [(The Arabidopsis Information Resource)](https://www.arabidopsis.org/) and [ARAPORT](https://apps.araport.org/thalemine/begin.do) contain the files necessary to do this.
