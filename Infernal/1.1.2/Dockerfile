FROM ubuntu:18.04
MAINTAINER Upendra Devisetty <upendra@cyverse.org>
LABEL version:1.1.2 Description:"This Dockerfile is for Infernal-1.1.2 tool"

RUN apt-get update && apt-get install -y wget make g++

RUN wget eddylab.org/infernal/infernal-1.1.2.tar.gz
RUN tar xf infernal-1.1.2.tar.gz
WORKDIR infernal-1.1.2
RUN ./configure && make && make install
WORKDIR easel
RUN make install
WORKDIR /
RUN wget ftp://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/Rfam.cm.gz
RUN gunzip Rfam.cm.gz

COPY Rfam.clanin /

RUN apt-get install -y bc

COPY infernal_wrapper.sh /usr/bin
RUN chmod +x /usr/bin/infernal_wrapper.sh

ENTRYPOINT ["infernal_wrapper.sh"]
