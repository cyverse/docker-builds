## Dockerfile
FROM ubuntu:18.04
MAINTAINER Amanda Cooksey	
LABEL Description="mugsy"

# Install all the updates and download dependencies
RUN apt-get update && \
    apt-get install -y \
    git \
    wget \
    bzip2

RUN echo 'export PATH=/opt/conda/bin:$PATH' > /etc/profile.d/conda.sh && \
    wget --quiet https://repo.continuum.io/miniconda/Miniconda3-4.7.12.1-Linux-x86_64.sh -O ~/miniconda.sh && \
    /bin/bash ~/miniconda.sh -b -p /opt/conda && \
    rm ~/miniconda.sh


# give write permissions to conda folder
RUN chmod 777 -R /opt/conda/

ENV PATH=$PATH:/opt/conda/bin

RUN conda config --add channels bioconda

RUN conda upgrade conda

RUN conda install -c conda-forge -c bioconda mugsy==1.2.3

ENV MUGSY_INSTALL=/opt/conda/bin

ENV PATH=$PATH:/opt/conda/bin/MUMmer3.20

WORKDIR /workdir
