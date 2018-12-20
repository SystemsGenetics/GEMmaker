# Base image.
FROM systemsgenetics/gemmaker:latest

USER root

# Install salmon verison 0.12.0
RUN wget https://github.com/COMBINE-lab/salmon/releases/download/v0.12.0/salmon-0.12.0_linux_x86_64.tar.gz \
  && tar xvzf salmon-0.12.0_linux_x86_64.tar.gz \
  && mv salmon-0.12.0_linux_x86_64 /usr/local/salmon-0.12.0 \
  && chmod 755 /usr/local/salmon-0.12.0/bin/salmon \
  && ln -sfn /usr/local/salmon-0.12.0/bin/salmon /usr/bin/salmon

USER gemdocker

