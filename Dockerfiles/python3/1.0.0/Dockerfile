# Base image.
FROM systemsgenetics/gemmaker:latest

USER root

RUN apt-get install -y python3-pip

RUN pip3 install --upgrade pip
RUN pip3 install numpy pandas

RUN git clone --single-branch -b issue_38_docker_integration https://github.com/SystemsGenetics/GEMmaker.git \
    && cp -r GEMmaker/scripts/*.py /usr/bin/

RUN  pip install -r GEMmaker/scripts/requirements.txt

USER gemdocker
