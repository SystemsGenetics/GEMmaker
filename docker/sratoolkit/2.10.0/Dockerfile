FROM systemsgenetics/aspera:3.8.1

# Install the SRA toolkit.
RUN wget -q https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/2.10.0/sratoolkit.2.10.0-ubuntu64.tar.gz \
  && tar -xf sratoolkit.2.10.0-ubuntu64.tar.gz \
  && mv sratoolkit.2.10.0-ubuntu64 /usr/local/sratoolkit.2.10.0

ENV PATH "$PATH:/usr/local/sratoolkit.2.10.0/bin"

# We need Python for the retrieve_sra.py script that comes with GEMmaker
RUN apt-get update && apt-get install -qq -y python3-pip \
  && pip3 install -q numpy pandas xmltodict

# The SRAToolkit needs a missing perl module
RUN apt-get install -qq -y libxml2-dev && cpan install -T XML::LibXML && cpan install -T URI
