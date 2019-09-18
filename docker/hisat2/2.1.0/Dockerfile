FROM gemmaker/base:1.1

# Install hisat2
RUN wget -q ftp://ftp.ccb.jhu.edu/pub/infphilo/hisat2/downloads/hisat2-2.1.0-Linux_x86_64.zip \
  && unzip -q hisat2-2.1.0-Linux_x86_64.zip \
  && mv hisat2-2.1.0 /usr/local

ENV PATH "$PATH:/usr/local/hisat2-2.1.0"
