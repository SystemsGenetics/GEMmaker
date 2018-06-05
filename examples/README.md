# Example Data

## Local Data Example Run
This directory contains both the reference material for the imaginary organism "Cool Organism" (CORG) and a dataset of 3 artificially made runs. This data set is very small, so it can be run on a local machine with ease.

The reference directory contains the made up reference genome file (**CORG.fna**), gtf file (**CORG.gtf**), and hisat files (**CORG.?/ht2**).

## Remote Data Example Run
This directory also contains a file titled "REMOTE_IDs.txt". This file is an example of a file that is in the proper format to download files from NCBI. To run this remote example, simpley change the "nextflow.config" file's "remote_list_path" parameter to be
```
'remote_list_path = "${PWD}/examples/REMOTE_IDs.txt"'
```
You can also change the "local_samples_path" to be "none" if you do not want the local files to run. This workflow can run both local and remote at the same time, however.

It should be noted that the output of the remote example is not as exciting as the local example, and is intended to be a demonstration of the download process works rather than the output process. The remote example uses the CORG reference files even though the remote data is *Arabidopsis thaliana*(Mouse Ear Crest). The final fpkm files will have values of 0 because of this, but the workflow will run as normal so that you can see how the remote process works. The reason that we have not included the *A. thaliana* reference files is to keep this git repository small.

It may be good practice for a new user of this workflow to generate the *A. thaliana* reference files that are associated with these small rna-seq datasets on their own. TAIR [(The Arabidopsis Information Resource)](https://www.arabidopsis.org/) and [ARAPORT](https://apps.araport.org/thalemine/begin.do) contain the files necessary to do this.
