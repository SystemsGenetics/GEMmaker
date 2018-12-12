FROM python:3.6


RUN apt-get update \
  && apt-get -y install git

RUN git clone --single-branch -b issue_38_docker_integration https://github.com/SystemsGenetics/GEMmaker.git \
    && cp -r GEMmaker/scripts/*.py /usr/bin/

RUN  pip install -r GEMmaker/scripts/requirements.txt

# Tini for signal processing and zombie killing.
ENV TINI_VERSION v0.18.0
ADD https://github.com/krallin/tini/releases/download/${TINI_VERSION}/tini /tini
RUN chmod +x /tini


# Create a non-root user that will run the called programs.
RUN groupadd -g 61000 gemdocker \
  && useradd -g 61000 --no-log-init --no-create-home --shell /bin/false -u 61000 gemdocker

# Change the active user.
USER gemdocker


# Define the command and parameters that will be executed when this
# container is first run.
ENTRYPOINT ["/tini", "--"]
