FROM ubuntu:16.04
MAINTAINER Upendra Devisetty <upendra@cyverse.org>
LABEL version="0.24" description="This Dockerfile is for process-reads-fasta"

RUN apt-get update && apt-get install -y wget python

# miR-PREFER
ENV miRPREFER_VERSION=0.24
RUN wget https://github.com/hangelwen/miR-PREFeR/archive/v${miRPREFER_VERSION}.tar.gz
RUN tar xvf v${miRPREFER_VERSION}.tar.gz

# Copy wrapper script
COPY process-reads-fasta_wrapper.sh /usr/bin
RUN chmod +x /usr/bin/process-reads-fasta_wrapper.sh

# ENTRYPOINT
ENTRYPOINT ["process-reads-fasta_wrapper.sh"]
