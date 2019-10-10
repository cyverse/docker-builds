FROM ubuntu:18.04
MAINTAINER Upendra Devisetty <upendra@cyverse.org>
LABEL Description "This Dockerfile is for SPRINT-0.1.8"

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
RUN conda install -y samtools==1.2 bwa==0.7.12

RUN wget -O- https://github.com/jumphone/SPRINT/archive/0.1.8.tar.gz | tar xzvf -

RUN chmod +x /SPRINT-0.1.8/bin/sprint 
RUN chmod +x /SPRINT-0.1.8/bin/sprint_from_bam 

COPY sprint_wrapper.sh /usr/bin
RUN chmod +x /usr/bin/sprint_wrapper.sh

ENV PATH /SPRINT-0.1.8/bin/:$PATH

ENTRYPOINT [ "sprint_wrapper.sh" ]
