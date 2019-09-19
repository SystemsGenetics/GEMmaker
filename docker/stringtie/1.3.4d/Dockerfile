FROM gemmaker/base:1.1

# Install stringtie.
RUN wget -q http://ccb.jhu.edu/software/stringtie/dl/stringtie-1.3.4d.Linux_x86_64.tar.gz \
  && tar -xf stringtie-1.3.4d.Linux_x86_64.tar.gz \
  && mv stringtie-1.3.4d.Linux_x86_64 /usr/local/stringtie-1.3.4d

# Add the prepDE.py script used for creating raw counts
RUN apt-get install -y python2.7 python-pip

RUN wget https://ccb.jhu.edu/software/stringtie/dl/prepDE.py \
  && chmod 755 prepDE.py \
  && mv prepDE.py /usr/local/stringtie-1.3.4d

ENV PATH "$PATH:/usr/local/stringtie-1.3.4d"
