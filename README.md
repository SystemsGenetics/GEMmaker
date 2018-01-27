# sra2gev

The sra2gev project is a [NextFlow](https://www.nextflow.io/) that downloads a set of samples from the NCBI Short Read Archive (SRA) and generates a file containing FPKM values for all genes in a genome annotation set. In other words, a Gene Expression Vector (GEV) for each sample is created. This workflow combines the [sratoolkit](https://www.ncbi.nlm.nih.gov/books/NBK158900/), [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic), [hisat2](https://ccb.jhu.edu/software/hisat2/index.shtml), [Samtools](http://www.htslib.org/), and [StringTie](http://www.ccb.jhu.edu/software/stringtie/).  The workflow expects the Lua-based [Lmod](https://lmod.readthedocs.io/en/latest/) software module system is installed with each of software described above availalbe via the module system.  The sra2gev workflow is setup to work with Illumina RNA-seq datasets housed in the SRA.

For testing purpose, or for execution of a small dataset (or large dataset if sufficient storage is available), a Docker image is available that contains all of the necessary software components: https://github.com/SystemsGenetics/sra2gev-docker

## Prepare the Workflow
Clone this workflow project into a working directory.  As with all NextFlow workflows, you can configure the behavior of the workflow by creating a **nextflow.config** file.  The sra2gev workflow provides an example file you can use to get started.  Rename this example to get started:

```bash
cp nextflow.config.example nextflow.config
```
To execute the sra2gev with the example dataset you must ensure that the **trimmomatic.clip_path** option is set to the full path where the Trimmomatic clipping files are housed.  Replace the text **<ILLUMINACLIP_PATH>** placeholder text. The example config file also has an example profile for running this workflow on a SLURM cluster. To use the SLURM profile you must, at a minimum, change the **<QUEUE_NAME>** placeholder text to be the name of the queue used for submission on your cluster.  If you require additional settings you can adjust the profile as per the [NextFlow configuration documentation](https://www.nextflow.io/docs/latest/config.html#config-profiles).  

## Execute the Workflow
To execute the workflow on a local machine use this command

```bash
nextflow run sra2gev.nf -profile standard
```

To resume a workflow in the event of a failure:
```bash
nextflow run sra2gev.nf -profile standard resume
```



