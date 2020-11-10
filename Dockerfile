FROM nfcore/base:1.11
LABEL authors="Simon Heumos, Michael Heuer, Erik Garrison, Andrea Guarracino" \
      description="Docker image containing all software requirements for the nf-core/pangenome pipeline"

# Install the conda environment
COPY environment.yml /
RUN conda env create --quiet -f /environment.yml && conda clean -a

# Copied from https://github.com/pangenome/pggb/blob/master/Dockerfile
RUN apt-get update \
    && apt-get install -y \
                       git \
                       bash \
                       cmake \
                       make \
                       g++ \
                       python3-dev \
                       libatomic-ops-dev

RUN cd ../../
RUN git clone --recursive https://github.com/ekg/edyeet
RUN apt-get install -y \
                        autoconf \
                        libgsl-dev \
                        zlib1g-dev
RUN cd edyeet \
    && git pull \
    && git checkout 1a172f6 \
    && bash bootstrap.sh \
    && bash configure \
    && make \
    && cp edyeet /usr/local/bin/edyeet

RUN cd ../
RUN git clone --recursive https://github.com/ekg/wfmash
RUN cd wfmash \
    && git pull \
    && git checkout 4bb309a \
    && bash bootstrap.sh \
    && bash configure \
    && make \
    && cp wfmash /usr/local/bin/wfmash

RUN cd ../
RUN git clone --recursive https://github.com/ekg/seqwish
RUN apt-get install -y \
                        build-essential
RUN cd seqwish \
    && git pull \
    && git checkout 9bbfa70 \
    && cmake -H. -Bbuild && cmake --build build -- -j $(nproc) \
    && cp bin/seqwish /usr/local/bin/seqwish

RUN cd ../
RUN git clone --recursive https://github.com/ekg/smoothxg
RUN cd smoothxg \
    && git pull \
    && git submodule update \
    && git checkout ffc50b7 \
    && sed -i 's/-march=native/-march=haswell/g' deps/abPOA/CMakeLists.txt \
    && cmake -H. -Bbuild && cmake --build build -- -j $(nproc) \
    && cp bin/smoothxg /usr/local/bin/smoothxg \
    && cp deps/odgi/bin/odgi /usr/local/bin/odgi

# Add conda installation dir to PATH (instead of doing 'conda activate')
ENV PATH /opt/conda/envs/nf-core-pangenome-1.0dev/bin:/usr/local/bin:$PATH

# Dump the details of the conda-installed packages to a file for posterity
RUN conda env export --name nf-core-pangenome-1.0dev > nf-core-pangenome-1.0dev.yml

# Instruct R processes to use these empty files instead of clashing with a local version
RUN touch .Rprofile
RUN touch .Renviron
