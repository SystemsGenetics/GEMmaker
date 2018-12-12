# Install FastQC
RUN wget https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.7.zip \
  && unzip fastqc_v0.11.7.zip \
  && mv FastQC /usr/local/FastQC-0.11.7 \
  && chmod 755 /usr/local/FastQC-0.11.7/fastqc \
  && ln -sfn /usr/local/FastQC-0.11.7/fastqc /usr/bin/fastqc
