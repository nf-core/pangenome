FROM ghcr.io/pangenome/pggb:20230113201558a9a04c
LABEL authors="Simon Heumos, Michael Heuer, Lukas Heumos, Erik Garrison, Andrea Guarracino" \
      description="Docker image containing all software requirements for the nf-core/pangenome pipeline"

# Install procps so that Nextflow can poll CPU usage and
# deep clean the apt cache to reduce image/layer size
RUN apt-get update \
    && apt-get install -y wget \
                          procps \
    && apt-get clean -y && rm -rf /var/lib/apt/lists/*

COPY bin/split_approx_mappings_in_chunks.py /
COPY bin/paf2net.py /
COPY bin/net2communities.py /

# Install miniconda
RUN wget \
    https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh \
    && mkdir /root/.conda \
    && bash Miniconda3-latest-Linux-x86_64.sh -b \
    && rm -f Miniconda3-latest-Linux-x86_64.sh

# Add conda to the path
ENV PATH /root/miniconda3/bin:$PATH

# Install the conda environment
COPY environment.yml /
RUN conda env create --quiet -f /environment.yml && conda clean -a

# Add conda installation dir to PATH (instead of doing 'conda activate')
ENV PATH /opt/conda/envs/nf-core-pangenome-1.0dev/bin:$PATH

# Set path for all users
RUN echo "export PATH=$PATH" > /etc/profile

# Dump the details of the conda-installed packages to a file for posterity
RUN conda env export --name nf-core-pangenome-1.0dev > nf-core-pangenome-1.0dev.yml

# this should fix https://nextflow.slack.com/archives/C02T98A23U7/p1671183108988419
SHELL ["/bin/bash", "-c"]
