FROM systemsgenetics/gemmaker:1.0

RUN apt-get install -qq -y autoconf automake make gcc perl zlib1g-dev libbz2-dev liblzma-dev libcurl4-gnutls-dev libssl-dev libncurses5-dev bzip2

# Install samtools
RUN wget -q https://github.com/samtools/samtools/releases/download/1.9/samtools-1.9.tar.bz2 \
  && bunzip2 -q samtools-1.9.tar.bz2 \
  && tar -xf samtools-1.9.tar \
  && cd samtools-1.9 \
  && ./configure --prefix=/usr/local/samtools-1.9 --without-curses \
  && make -s install

# Add those tools to the path variable.
ENV PATH="/usr/local/samtools-1.9/bin:${PATH}"
