FROM gemmaker/base:1.0

RUN apt-get install -qq -y default-jre perl

# Install FastQC
RUN wget -q https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.7.zip \
  && unzip -q fastqc_v0.11.7.zip \
  && mv FastQC /usr/local/FastQC-0.11.7 \
  && chmod 755 /usr/local/FastQC-0.11.7/fastqc \
  && ln -sfn /usr/local/FastQC-0.11.7/fastqc /usr/bin/fastqc
