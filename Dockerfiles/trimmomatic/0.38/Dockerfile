FROM systemsgenetics/gemmaker:latest

RUN apt-get update \
  && apt-get -y install default-jre procps unzip wget

# Install trimmomatic
RUN wget http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.38.zip \
  && unzip Trimmomatic-0.38.zip \
  && mv Trimmomatic-0.38 /usr/local

ENV PATH "$PATH:/usr/local/Trimmomatic-0.38/"
ENV CLASSPATH "/usr/local/Trimmomatic-0.38/trimmomatic-0.38.jar"
ENV ILLUMINACLIP_PATH "/usr/local/Trimmomatic-0.38/adapters"
