FROM systemsgenetics/gemmaker:latest

# Install kallisto verison 0.45.0
RUN wget https://github.com/pachterlab/kallisto/releases/download/v0.45.0/kallisto_linux-v0.45.0.tar.gz \
  && tar xvzf kallisto_linux-v0.45.0.tar.gz  \
  && mv kallisto_linux-v0.45.0 /usr/local/kallisto-0.45.0 \
  && chmod 755 /usr/local/kallisto-0.45.0/kallisto \
  && ln -sfn /usr/local/kallisto-0.45.0/kallisto /usr/bin/kallisto
