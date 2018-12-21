FROM systemsgenetics/gemmaker:latest

# Install hisat2
RUN wget ftp://ftp.ccb.jhu.edu/pub/infphilo/hisat2/downloads/hisat2-2.1.0-Linux_x86_64.zip \
  && unzip hisat2-2.1.0-Linux_x86_64.zip \
  && mv hisat2-2.1.0 /usr/local

ENV PATH "$PATH:/usr/local/hisat2-2.1.0"
