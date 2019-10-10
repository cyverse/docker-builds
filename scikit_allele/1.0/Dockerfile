FROM debian:8.5
MAINTAINER Upendra Kumar Devisetty <upendra@cyverse.org>
LABEL Description "This Dockerfile is for scikit-allel"

ENV LANG=C.UTF-8 LC_ALL=C.UTF-8

RUN apt-get update --fix-missing && apt-get install -y wget bzip2 ca-certificates \
    libglib2.0-0 libxext6 libsm6 libxrender1 \
    git mercurial subversion

RUN echo 'export PATH=/opt/conda/bin:$PATH' > /etc/profile.d/conda.sh && \
    wget --quiet https://repo.continuum.io/miniconda/Miniconda2-4.3.14-Linux-x86_64.sh -O ~/miniconda.sh && \
    /bin/bash ~/miniconda.sh -b -p /opt/conda && \
    rm ~/miniconda.sh

ENV PATH /opt/conda/bin:$PATH

RUN conda config --add channels conda-forge 
RUN conda config --add channels defaults 
RUN conda install scikit-allel

ADD extract_genotype.py /usr/bin
RUN chmod +x /usr/bin/extract_genotype.py

ENTRYPOINT ["python", "/usr/bin/extract_genotype.py"]
