FROM systemsgenetics/gemmaker:latest

RUN apt-get install -qq -y python3-pip \
  && pip3 install -q numpy pandas

RUN git clone --single-branch -b issue_38_docker_integration https://github.com/SystemsGenetics/GEMmaker.git \
  && pip3 install -q -r GEMmaker/bin/requirements.txt \
  && cp GEMmaker/bin/*.py /usr/bin/
