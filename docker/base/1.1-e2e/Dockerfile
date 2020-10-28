FROM ubuntu:18.04



# install apt dependencies
RUN apt-get update -qq \
  && apt-get install -qq -y \
    curl git unzip wget \
    default-jre perl \
    python3-pip \
    autoconf automake make gcc perl zlib1g-dev libbz2-dev liblzma-dev libcurl4-gnutls-dev libssl-dev libncurses5-dev bzip2 \
    python2.7 python-pip \
    default-jre procps



# install aspera 3.8.1
RUN groupadd -g 61000 gemmaker \
  && useradd -g 61000 --no-log-init --create-home --shell /bin/bash -u 61000 gemmaker

USER gemmaker
WORKDIR /home/gemmaker

RUN wget -q https://download.asperasoft.com/download/sw/connect/3.8.1/ibm-aspera-connect-3.8.1.161274-linux-g2.12-64.tar.gz \
  && tar -xf ibm-aspera-connect-3.8.1.161274-linux-g2.12-64.tar.gz \
  && ./ibm-aspera-connect-3.8.1.161274-linux-g2.12-64.sh \
  && rm ibm-aspera-connect-3.8.1.161274-linux-g2.12-64.sh

USER root
WORKDIR /root

RUN mv /home/gemmaker/.aspera /opt/aspera

ENV PATH "$PATH:/opt/aspera/connect/bin"
ENV ASPERA_KEY "/opt/aspera/connect/etc/asperaweb_id_dsa.openssh"



# install FastQC 0.11.7
RUN wget -q https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.7.zip \
  && unzip -q fastqc_v0.11.7.zip \
  && mv FastQC /usr/local/FastQC-0.11.7 \
  && chmod 755 /usr/local/FastQC-0.11.7/fastqc \
  && ln -sfn /usr/local/FastQC-0.11.7/fastqc /usr/bin/fastqc



# install hisat2 2.1.0
RUN wget -q ftp://ftp.ccb.jhu.edu/pub/infphilo/hisat2/downloads/hisat2-2.1.0-Linux_x86_64.zip \
  && unzip -q hisat2-2.1.0-Linux_x86_64.zip \
  && mv hisat2-2.1.0 /usr/local

ENV PATH "$PATH:/usr/local/hisat2-2.1.0"



# install kallisto 0.45.0
RUN wget -q https://github.com/pachterlab/kallisto/releases/download/v0.45.0/kallisto_linux-v0.45.0.tar.gz \
  && tar -xf kallisto_linux-v0.45.0.tar.gz  \
  && mv kallisto_linux-v0.45.0 /usr/local/kallisto-0.45.0 \
  && chmod 755 /usr/local/kallisto-0.45.0/kallisto \
  && ln -sfn /usr/local/kallisto-0.45.0/kallisto /usr/bin/kallisto



# Install MultiQC
ENV LC_ALL=C.UTF-8
ENV LANG=C.UTF-8

RUN pip3 install -q https://github.com/ewels/MultiQC/archive/v1.8.zip



# install python 3.7 modules
RUN pip3 install -q decorator>=4.3.0 networkx numpy pandas requests xmltodict



# install salmon 0.12.0
RUN wget -q https://github.com/COMBINE-lab/salmon/releases/download/v0.12.0/salmon-0.12.0_linux_x86_64.tar.gz \
  && tar -xf salmon-0.12.0_linux_x86_64.tar.gz \
  && mv salmon-0.12.0_linux_x86_64 /usr/local/salmon-0.12.0 \
  && chmod 755 /usr/local/salmon-0.12.0/bin/salmon \
  && ln -sfn /usr/local/salmon-0.12.0/bin/salmon /usr/bin/salmon



# install samtools 1.9
RUN wget -q https://github.com/samtools/samtools/releases/download/1.9/samtools-1.9.tar.bz2 \
  && bunzip2 -q samtools-1.9.tar.bz2 \
  && tar -xf samtools-1.9.tar \
  && cd samtools-1.9 \
  && ./configure --prefix=/usr/local/samtools-1.9 --without-curses \
  && make -s install

ENV PATH="$PATH:/usr/local/samtools-1.9/bin"



# install sratoolkit 2.10.0
RUN wget -q https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/2.10.0/sratoolkit.2.10.0-ubuntu64.tar.gz \
  && tar -xf sratoolkit.2.10.0-ubuntu64.tar.gz \
  && mv sratoolkit.2.10.0-ubuntu64 /usr/local/sratoolkit.2.10.0

ENV PATH "$PATH:/usr/local/sratoolkit.2.10.0/bin"

# The SRAToolkit needs a missing perl module
RUN cpan install -T URI XML::LibXML



# install stringtie 1.3.4d
RUN wget -q http://ccb.jhu.edu/software/stringtie/dl/stringtie-1.3.4d.Linux_x86_64.tar.gz \
  && tar -xf stringtie-1.3.4d.Linux_x86_64.tar.gz \
  && mv stringtie-1.3.4d.Linux_x86_64 /usr/local/stringtie-1.3.4d



# install the prepDE.py script used for creating raw counts
RUN wget -q https://ccb.jhu.edu/software/stringtie/dl/prepDE.py \
  && chmod 755 prepDE.py \
  && mv prepDE.py /usr/local/stringtie-1.3.4d

ENV PATH "$PATH:/usr/local/stringtie-1.3.4d"



# install trimmomatic 0.38
RUN wget -q http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.38.zip \
  && unzip -q Trimmomatic-0.38.zip \
  && mv Trimmomatic-0.38 /usr/local

ENV PATH "$PATH:/usr/local/Trimmomatic-0.38/"
ENV CLASSPATH "/usr/local/Trimmomatic-0.38/trimmomatic-0.38.jar"
ENV ILLUMINACLIP_PATH "/usr/local/Trimmomatic-0.38/adapters"
