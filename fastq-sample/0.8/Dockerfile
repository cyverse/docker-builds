FROM ubuntu:18.04
MAINTAINER Upendra Devisetty <upendra@cyverse.org>

# Install dependencies
RUN apt-get update && apt-get install -y  wget gcc make libpcre3-dev libz-dev 

# Download source code
RUN wget http://homes.cs.washington.edu/~dcjones/fastq-tools/fastq-tools-0.8.tar.gz
RUN tar xvf fastq-tools-0.8.tar.gz
WORKDIR fastq-tools-0.8

# Compile the source code
RUN ./configure
RUN  make install fastq==0.8

# Set the Entrypoint
ENTRYPOINT [ "fastq-sample" ]
