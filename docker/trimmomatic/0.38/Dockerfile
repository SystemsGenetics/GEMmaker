FROM gemmaker/base:1.1

RUN apt-get update -qq \
  && apt-get -qq -y install default-jre procps

# Install trimmomatic
RUN wget -q http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.38.zip \
  && unzip -q Trimmomatic-0.38.zip \
  && mv Trimmomatic-0.38 /usr/local

ENV PATH "$PATH:/usr/local/Trimmomatic-0.38/"
ENV CLASSPATH "/usr/local/Trimmomatic-0.38/trimmomatic-0.38.jar"
ENV ILLUMINACLIP_PATH "/usr/local/Trimmomatic-0.38/adapters"
