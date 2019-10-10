FROM debian:8.5
MAINTAINER Upendra Kumar Devisetty <upendra@cyverse.org>
LABEL Description "This Dockerfile is for Fastqc_screen_v0.11.1"

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
RUN conda config --add channels r 
RUN conda config --add channels bioconda

RUN conda install fastq-screen==0.11.1 perl-gd perl-gdgraph perl-gdgraph-histogram bwa==0.7.15

ADD fastq_screen_wrapper.sh /usr/bin/
RUN chmod +x /usr/bin/fastq_screen_wrapper.sh

ENTRYPOINT ["fastq_screen_wrapper.sh"]
CMD [ "-h" ]