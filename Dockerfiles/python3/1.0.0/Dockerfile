FROM systemsgenetics/gemmaker:latest

RUN apt-get install -y python3-pip
  && pip3 install --upgrade pip
  && pip3 install numpy pandas

RUN git clone --single-branch -b issue_38_docker_integration https://github.com/SystemsGenetics/GEMmaker.git \
  && pip install -r GEMmaker/bin/requirements.txt
  && cp -r GEMmaker/bin/*.py /usr/bin/
