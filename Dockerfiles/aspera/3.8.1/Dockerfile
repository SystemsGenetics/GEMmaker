FROM systemsgenetics/gemmaker:latest

RUN groupadd -g 61000 gemdocker \
  && useradd -g 61000 --no-log-init --create-home --shell /bin/bash -u 61000 gemdocker
USER gemdocker
WORKDIR /home/gemdocker

RUN wget https://download.asperasoft.com/download/sw/connect/3.8.1/ibm-aspera-connect-3.8.1.161274-linux-g2.12-64.tar.gz \
  && tar -zxvf ibm-aspera-connect-3.8.1.161274-linux-g2.12-64.tar.gz \
  && ./ibm-aspera-connect-3.8.1.161274-linux-g2.12-64.sh \
  && rm ibm-aspera-connect-3.8.1.161274-linux-g2.12-64.sh

USER root
WORKDIR /root

# Make sure the ascp command is in the path and to support future
# updates of aspera we'll add a new variable ASPERA_KEY that can be
# used in the GEMmaker bash code.
ENV PATH "$PATH:/home/gemdocker/.aspera/connect/bin"
