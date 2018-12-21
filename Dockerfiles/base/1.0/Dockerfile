FROM ubuntu:16.04
MAINTAINER Tyler Biggs <tyler.biggs@wsu.edu>

LABEL base_image="ubuntu:16.04"
LABEL software="GEM Maker Base Image"
LABEL software.version="1.0.0"
LABEL about.summary="Base image for GEM Maker"
LABEL about.home=""
LABEL about.documentation=""
LABEL about.license_file=""

# Set noninteractive mode for apt-get.
ENV DEBIAN_FRONTEND noninteractive

# Update and install basic packages.
RUN apt-get update -qq \
  && apt-get install -qq -y curl git unzip wget 
