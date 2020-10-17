GEMMAKER_VERSION="1.1"

all:
	docker build -t gemmaker/base:${GEMMAKER_VERSION}             base/${GEMMAKER_VERSION}/
	docker build -t gemmaker/base:${GEMMAKER_VERSION}-merged      base/${GEMMAKER_VERSION}-merged/
	docker build -t gemmaker/aspera:3.8.1-${GEMMAKER_VERSION}     aspera/3.8.1-${GEMMAKER_VERSION}/
	docker build -t gemmaker/fastqc:0.11.7-${GEMMAKER_VERSION}    fastqc/0.11.7-${GEMMAKER_VERSION}/
	docker build -t gemmaker/hisat2:2.1.0-${GEMMAKER_VERSION}     hisat2/2.1.0-${GEMMAKER_VERSION}/
	docker build -t gemmaker/kallisto:0.45.0-${GEMMAKER_VERSION}  kallisto/0.45.0-${GEMMAKER_VERSION}/
	docker build -t gemmaker/python:3.5.1-${GEMMAKER_VERSION}     python/3.5.1-${GEMMAKER_VERSION}/
	docker build -t gemmaker/multiqc:1.7-${GEMMAKER_VERSION}      multiqc/1.7-${GEMMAKER_VERSION}/
	docker build -t gemmaker/salmon:0.12.0-${GEMMAKER_VERSION}    salmon/0.12.0-${GEMMAKER_VERSION}/
	docker build -t gemmaker/samtools:1.9-${GEMMAKER_VERSION}     samtools/1.9-${GEMMAKER_VERSION}/
	docker build -t gemmaker/sratoolkit:2.9.2-${GEMMAKER_VERSION} sratoolkit/2.9.2-${GEMMAKER_VERSION}/
	docker build -t gemmaker/sratoolkit:2.9.6-${GEMMAKER_VERSION} sratoolkit/2.9.6-${GEMMAKER_VERSION}/
	docker build -t gemmaker/stringtie:1.3.4d-${GEMMAKER_VERSION} stringtie/1.3.4d-${GEMMAKER_VERSION}/
	docker build -t gemmaker/trimmomatic:0.38-${GEMMAKER_VERSION} trimmomatic/0.38-${GEMMAKER_VERSION}/

push:
	docker push gemmaker/base:${GEMMAKER_VERSION}
	docker push gemmaker/base:${GEMMAKER_VERSION}-merged
	docker push gemmaker/aspera:3.8.1-${GEMMAKER_VERSION}
	docker push gemmaker/fastqc:0.11.7-${GEMMAKER_VERSION}
	docker push gemmaker/hisat2:2.1.0-${GEMMAKER_VERSION}
	docker push gemmaker/kallisto:0.45.0-${GEMMAKER_VERSION}
	docker push gemmaker/python:3.5.1-${GEMMAKER_VERSION}
	docker push gemmaker/multiqc:1.7-${GEMMAKER_VERSION}
	docker push gemmaker/salmon:0.12.0-${GEMMAKER_VERSION}
	docker push gemmaker/samtools:1.9-${GEMMAKER_VERSION}
	docker push gemmaker/sratoolkit:2.9.2-${GEMMAKER_VERSION}
	docker push gemmaker/sratoolkit:2.9.6-${GEMMAKER_VERSION}
	docker push gemmaker/stringtie:1.3.4d-${GEMMAKER_VERSION}
	docker push gemmaker/trimmomatic:0.38-${GEMMAKER_VERSION}

pull:
	docker pull gemmaker/base:${GEMMAKER_VERSION}
	docker pull gemmaker/base:${GEMMAKER_VERSION}-merged
	docker pull gemmaker/aspera:3.8.1-${GEMMAKER_VERSION}
	docker pull gemmaker/fastqc:0.11.7-${GEMMAKER_VERSION}
	docker pull gemmaker/hisat2:2.1.0-${GEMMAKER_VERSION}
	docker pull gemmaker/kallisto:0.45.0-${GEMMAKER_VERSION}
	docker pull gemmaker/python:3.5.1-${GEMMAKER_VERSION}
	docker pull gemmaker/multiqc:1.7-${GEMMAKER_VERSION}
	docker pull gemmaker/salmon:0.12.0-${GEMMAKER_VERSION}
	docker pull gemmaker/samtools:1.9-${GEMMAKER_VERSION}
	docker pull gemmaker/sratoolkit:2.9.2-${GEMMAKER_VERSION}
	docker pull gemmaker/sratoolkit:2.9.6-${GEMMAKER_VERSION}
	docker pull gemmaker/stringtie:1.3.4d-${GEMMAKER_VERSION}
	docker pull gemmaker/trimmomatic:0.38-${GEMMAKER_VERSION}

clean:
	docker image rm -f gemmaker/base:${GEMMAKER_VERSION}
	docker image rm -f gemmaker/base:${GEMMAKER_VERSION}-merged
	docker image rm -f gemmaker/aspera:3.8.1-${GEMMAKER_VERSION}
	docker image rm -f gemmaker/fastqc:0.11.7-${GEMMAKER_VERSION}
	docker image rm -f gemmaker/hisat2:2.1.0-${GEMMAKER_VERSION}
	docker image rm -f gemmaker/kallisto:0.45.0-${GEMMAKER_VERSION}
	docker image rm -f gemmaker/python:3.5.1-${GEMMAKER_VERSION}
	docker image rm -f gemmaker/multiqc:1.7-${GEMMAKER_VERSION}
	docker image rm -f gemmaker/salmon:0.12.0-${GEMMAKER_VERSION}
	docker image rm -f gemmaker/samtools:1.9-${GEMMAKER_VERSION}
	docker image rm -f gemmaker/sratoolkit:2.9.2-${GEMMAKER_VERSION}
	docker image rm -f gemmaker/sratoolkit:2.9.6-${GEMMAKER_VERSION}
	docker image rm -f gemmaker/stringtie:1.3.4d-${GEMMAKER_VERSION}
	docker image rm -f gemmaker/trimmomatic:0.38-${GEMMAKER_VERSION}
