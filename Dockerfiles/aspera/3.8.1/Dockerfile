FROM systemsgenetics/sratoolkit:2.9.6

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

# Make sure the ascp command is in the path and to support future
# updates of aspera we'll add a new variable ASPERA_KEY that can be
# used in the GEMmaker bash code.
ENV PATH "$PATH:/opt/aspera/connect/bin"
ENV ASPERA_KEY "/opt/aspera/connect/etc/asperaweb_id_dsa.openssh"
