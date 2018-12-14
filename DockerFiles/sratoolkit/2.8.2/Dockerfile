# Base image.
FROM systemsgenetics/gemmaker:latest

USER root

ENV SRA_VERSION 2.8.2

# Install the SRA toolkit.
RUN wget -q https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/$SRA_VERSION/sratoolkit.$SRA_VERSION-ubuntu64.tar.gz \
    && tar -zxvf sratoolkit.$SRA_VERSION-ubuntu64.tar.gz \
    && mv sratoolkit.$SRA_VERSION-ubuntu64 /usr/local/sratoolkit.$SRA_VERSION

USER gemdocker
