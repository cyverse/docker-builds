## Dockerfile
FROM ubuntu:18.04
MAINTAINER Amanda Cooksey	
LABEL Description="vcf2phylip"

# Install all the updates and download dependencies
RUN apt-get update && \
    apt-get install -y \
    git \
    wget \
    bzip2 \
    python3 \
    nano

WORKDIR /usr/local/bin/

RUN git clone https://github.com/edgardomortiz/vcf2phylip.git 

WORKDIR /workdir

ENV PATH=$PATH:/usr/local/bin/vcf2phylip



