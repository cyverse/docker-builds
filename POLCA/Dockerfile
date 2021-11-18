## Dockerfile
FROM ubuntu:20.04
MAINTAINER Amanda Cooksey	
LABEL Description="POLCA"

# Install all the updates and download dependencies
RUN apt-get update && \
    apt-get install -y \
    git \
    wget \
    bzip2 

ENV PATH="/root/miniconda3/bin:${PATH}"
ARG PATH="/root/miniconda3/bin:${PATH}"
RUN wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh \
    && mkdir /root/.conda \
    && bash Miniconda3-latest-Linux-x86_64.sh -b \
    && rm -f Miniconda3-latest-Linux-x86_64.sh


RUN conda config --add channels bioconda

RUN conda upgrade conda

RUN conda install -c conda-forge -c bioconda masurca==3.4.2 
RUN conda install -c conda-forge -c bioconda bwa==0.7.17

ENTRYPOINT ["polca.sh"]
