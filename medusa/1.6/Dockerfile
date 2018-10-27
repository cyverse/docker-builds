FROM ubuntu:16.04
MAINTAINER Upendra Devisetty <upendra@cyverse.org>
LABEL version="1.6" description="This Dockerfile is for medusa"

RUN apt-get update && \
    apt-get install -y default-jre && \ 
    apt-get install -y python-pip && \
    pip install networkx && \
    pip install biopython

RUN apt-get install -y wget

RUN wget https://github.com/mummer4/mummer/releases/download/v4.0.0beta2/mummer-4.0.0beta2.tar.gz
RUN tar xvf mummer-4.0.0beta2.tar.gz
WORKDIR mummer-4.0.0beta2
RUN ./configure --prefix=/mummer-4.0.0beta2 && make && make install
ENV PATH=/mummer-4.0.0beta2/bin/:$PATH

WORKDIR /										 
RUN wget https://github.com/combogenomics/medusa/releases/download/medusa-1.6/medusa.tar.gz
RUN tar xvf medusa.tar.gz

ENTRYPOINT ["java", "-jar", "/medusa/medusa.jar"]
