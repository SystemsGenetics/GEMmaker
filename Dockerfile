FROM nfcore/base:1.13.3
LABEL authors="John Hadish, Tyler Biggs, Ben Shealy, Connor Wytko, Sai Prudhvi Oruganti, F. Alex Feltus, & Stephen Ficklin" \
      description="Docker image containing all software requirements for the nf-core/gemmaker pipeline"

# Install the conda environment
COPY environment.yml /
RUN conda env create --quiet -f /environment.yml && conda clean -a

# Add conda installation dir to PATH (instead of doing 'conda activate')
ENV PATH /opt/conda/envs/nf-core-gemmaker-2.0dev/bin:$PATH

# Dump the details of the installed packages to a file for posterity
RUN conda env export --name nf-core-gemmaker-2.0dev > nf-core-gemmaker-2.0dev.yml
