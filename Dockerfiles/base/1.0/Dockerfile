# Base image.
FROM ubuntu:16.04


################## METADATA ######################
LABEL base_image="ubuntu:16.04"
LABEL software="GEM Maker Base Image"
LABEL software.version="1.0.0"
LABEL about.summary="Base image for GEM Maker"
LABEL about.home=""
LABEL about.documentation=""
LABEL about.license_file=""


################## MAINTAINER ######################
MAINTAINER Tyler Biggs <tyler.biggs@wsu.edu>


# Set noninteractive mode for apt-get.
ENV DEBIAN_FRONTEND noninteractive


# Update and install basic packages.
RUN apt-get update && apt-get upgrade -y wget curl unzip \
  && apt-get install -y git


# Create a non-root user that will run the called programs.
RUN groupadd -g 61000 gemdocker \
  && useradd -g 61000 --no-log-init --no-create-home --shell /bin/bash -u 61000 gemdocker

# Change the active user.
USER gemdocker
