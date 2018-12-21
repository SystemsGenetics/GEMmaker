FROM systemsgenetics/gemmaker:latest

# Install the SRA toolkit.
RUN wget -q https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/2.9.2/sratoolkit.2.9.2-ubuntu64.tar.gz \
  && tar -zxvf sratoolkit.2.9.2-ubuntu64.tar.gz \
  && mv sratoolkit.2.9.2-ubuntu64 /usr/local/sratoolkit.2.9.2

ENV PATH "$PATH:/usr/local/sratoolkit.2.9.2/bin"
