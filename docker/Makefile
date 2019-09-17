all:
	docker build -t systemsgenetics/gemmaker:1.0        gemmaker/1.0/
	docker build -t systemsgenetics/gemmaker:1.0-merged gemmaker/1.0-merged/
	docker build -t systemsgenetics/aspera:3.8.1        aspera/3.8.1/
	docker build -t systemsgenetics/fastqc:0.11.7       fastqc/0.11.7/
	docker build -t systemsgenetics/hisat2:2.1.0        hisat2/2.1.0/
	docker build -t systemsgenetics/kallisto:0.45.0     kallisto/0.45.0/
	docker build -t systemsgenetics/python:3.5.1        python/3.5.1/
	docker build -t systemsgenetics/multiqc:1.7         multiqc/1.7
	docker build -t systemsgenetics/salmon:0.12.0       salmon/0.12.0/
	docker build -t systemsgenetics/samtools:1.9        samtools/1.9/
	docker build -t systemsgenetics/sratoolkit:2.9.2    sratoolkit/2.9.2/
	docker build -t systemsgenetics/sratoolkit:2.9.6    sratoolkit/2.9.6/
	docker build -t systemsgenetics/stringtie:1.3.4d    stringtie/1.3.4d/
	docker build -t systemsgenetics/trimmomatic:0.38    trimmomatic/0.38/

push:
	docker push systemsgenetics/gemmaker:1.0
	docker push systemsgenetics/gemmaker:1.0-merged
	docker push systemsgenetics/aspera:3.8.1
	docker push systemsgenetics/fastqc:0.11.7
	docker push systemsgenetics/hisat2:2.1.0
	docker push systemsgenetics/kallisto:0.45.0
	docker push systemsgenetics/python:3.5.1
	docker push systemsgenetics/multiqc:1.7
	docker push systemsgenetics/salmon:0.12.0
	docker push systemsgenetics/samtools:1.9
	docker push systemsgenetics/sratoolkit:2.9.2
	docker push systemsgenetics/sratoolkit:2.9.6
	docker push systemsgenetics/stringtie:1.3.4d
	docker push systemsgenetics/trimmomatic:0.38

pull:
	docker pull systemsgenetics/gemmaker:1.0
	docker pull systemsgenetics/gemmaker:1.0-merged
	docker pull systemsgenetics/aspera:3.8.1
	docker pull systemsgenetics/fastqc:0.11.7
	docker pull systemsgenetics/hisat2:2.1.0
	docker pull systemsgenetics/kallisto:0.45.0
	docker pull systemsgenetics/python:3.5.1
	docker pull systemsgenetics/multiqc:1.7
	docker pull systemsgenetics/salmon:0.12.0
	docker pull systemsgenetics/samtools:1.9
	docker pull systemsgenetics/sratoolkit:2.9.2
	docker pull systemsgenetics/sratoolkit:2.9.6
	docker pull systemsgenetics/stringtie:1.3.4d
	docker pull systemsgenetics/trimmomatic:0.38

clean:
	docker image rm -f systemsgenetics/gemmaker:1.0
	docker image rm -f systemsgenetics/gemmaker:1.0-merged
	docker image rm -f systemsgenetics/aspera:3.8.1
	docker image rm -f systemsgenetics/fastqc:0.11.7
	docker image rm -f systemsgenetics/hisat2:2.1.0
	docker image rm -f systemsgenetics/kallisto:0.45.0
	docker image rm -f systemsgenetics/python:3.5.1
	docker image rm -f systemsgenetics/multiqc:1.7
	docker image rm -f systemsgenetics/salmon:0.12.0
	docker image rm -f systemsgenetics/samtools:1.9
	docker image rm -f systemsgenetics/sratoolkit:2.9.2
	docker image rm -f systemsgenetics/sratoolkit:2.9.6
	docker image rm -f systemsgenetics/stringtie:1.3.4d
	docker image rm -f systemsgenetics/trimmomatic:0.38
