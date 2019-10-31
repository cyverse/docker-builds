FROM ubuntu:18.04
MAINTAINER Upendra Devisetty <upendra@cyverse.org>
LABEL Description "This Dockerfile is for Trimmomatic-0.38"

RUN apt-get update && apt-get install -y wget python2.7

ENV LANG=C.UTF-8 LC_ALL=C.UTF-8

RUN echo 'export PATH=/opt/conda/bin:$PATH' > /etc/profile.d/conda.sh && \
    wget --quiet https://repo.anaconda.com/miniconda/Miniconda2-latest-Linux-x86_64.sh -O ~/miniconda.sh && \
    /bin/bash ~/miniconda.sh -b -p /opt/conda && \
    rm ~/miniconda.sh

ENV PATH /opt/conda/bin:$PATH

RUN conda config --add channels conda-forge 
RUN conda config --add channels defaults 
RUN conda config --add channels r 
RUN conda config --add channels bioconda
RUN conda install -y trimmomatic==0.38 

COPY trimscript.wrapper_0.38.sh /usr/bin
RUN chmod +x /usr/bin/trimscript.wrapper_0.38.sh

ENTRYPOINT [ "trimscript.wrapper_0.38.sh" ]
