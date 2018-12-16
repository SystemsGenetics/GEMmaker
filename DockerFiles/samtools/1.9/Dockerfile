# Base image.
FROM systemsgenetics/gemmaker:latest

USER root

RUN apt-get install -y autoconf automake make gcc perl zlib1g-dev libbz2-dev liblzma-dev libcurl4-gnutls-dev libssl-dev libncurses5-dev bzip2

# Install samtools
RUN wget https://github.com/samtools/samtools/releases/download/1.9/samtools-1.9.tar.bz2 \
  && bunzip2 samtools-1.9.tar.bz2 \
  && tar -xvf samtools-1.9.tar \
  && cd samtools-1.9 \
  && ./configure --prefix=/usr/local/samtools-1.9 --without-curses \
  && make install
# Add those tools to the path variable.
ENV PATH="/usr/local/samtools-1.9/bin:${PATH}"


USER gemdocker
