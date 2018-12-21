FROM systemsgenetics/gemmaker:latest

# Install stringtie.
RUN wget -q http://ccb.jhu.edu/software/stringtie/dl/stringtie-1.3.4d.Linux_x86_64.tar.gz \
  && tar -xf stringtie-1.3.4d.Linux_x86_64.tar.gz \
  && mv stringtie-1.3.4d.Linux_x86_64 /usr/local/stringtie-1.3.4d

ENV PATH "$PATH:/usr/local/stringtie-1.3.4d"
