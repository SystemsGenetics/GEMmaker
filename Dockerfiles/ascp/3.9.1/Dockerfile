FROM systemsgenetics/gemmaker:1.0

# Version 3.9.1
# Platform: Linux x86_64
# Release Date: December 7, 2018
#
# md5: 743da65c7adbe48558a6890bb1bcb007
# sha1: edeaa6d70e5b25e5d092b047c8429fa0c5ae9bda
USER root

RUN wget https://download.asperasoft.com/download/sw/client/3.9.1/ibm-aspera-desktopclient-3.9.1.168302-linux-64.deb \
  && dpkg -i ibm-aspera-desktopclient-3.9.1.168302-linux-64.deb \
  && apt-get purge -y wget \
  && apt-get autoremove -y \
  && apt-get clean && \
  rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

USER gemdocker
