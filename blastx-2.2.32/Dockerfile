# Pull base image
FROM ubuntu:12.04

MAINTAINER Andre Mercer <amercer@iplantcollaborative.org>

# Run updates
RUN apt-get update
RUN apt-get upgrade -y
RUN apt-get install wget -y

# Pull tool
RUN wget -P /root ftp://ftp.ncbi.nlm.nih.gov/blast/executables/LATEST/ncbi-blast-2.2.31+-x64-linux.tar.gz
# Open tar
RUN tar zxf /root/ncbi-blast-2.2.31+-x64-linux.tar.gz -C /root

RUN mkdir /root/ncbi-blast-2.2.31+/db

ENV PATH "/root/ncbi-blast-2.2.31+/bin:$PATH"
ENV BLASTDB "/root/ncbi-blast-2.2.31+/db"

ENTRYPOINT ["blastx"]
