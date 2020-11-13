FROM ghcr.io/pangenome/pggb:latest
LABEL authors="Simon Heumos, Michael Heuer, Lukas Heumos, Erik Garrison, Andrea Guarracino" \
      description="Docker image containing all software requirements for the nf-core/pangenome pipeline"

# Install procps so that Nextflow can poll CPU usage and
# deep clean the apt cache to reduce image/layer size
RUN apt-get update \
    && apt-get install -y wget \
                          procps \
    && apt-get clean -y && rm -rf /var/lib/apt/lists/*

# Install Miniconda
RUN wget -O ~/miniconda.sh https://repo.continuum.io/miniconda/Miniconda3-py37_4.8.2-Linux-x86_64.sh \
 && chmod +x ~/miniconda.sh \
 && ~/miniconda.sh -b -p ~/miniconda \
 && rm ~/miniconda.sh
ENV PATH "$PATH:/root/miniconda/bin"
RUN echo "export PATH=$PATH" > /etc/profile
ENV CONDA_AUTO_UPDATE_CONDA=false

# Install the conda environment
COPY environment.yml /
RUN conda env create --quiet -f /environment.yml && conda clean -a

# Add conda installation dir to runtime PATH (instead of doing 'conda activate')
ENV PATH /opt/conda/envs/nf-core-pangenome-1.0dev/bin:$PATH

# Dump the details of the conda-installed packages to a file for posterity
RUN conda env export --name nf-core-pangenome-1.0dev > nf-core-pangenome-1.0dev.yml

# Instruct R processes to use these empty files instead of clashing with a local version
RUN touch .Rprofile
RUN touch .Renviron
