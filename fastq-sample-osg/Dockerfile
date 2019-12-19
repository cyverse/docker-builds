FROM ubuntu:xenial
MAINTAINER Upendra Devisetty <upendra@cyverse.org>

RUN mkdir /cvmfs /work

RUN apt-get update \
    && apt-get install -y lsb curl apt-transport-https python3 python-requests libfuse2 wget gcc make libpcre3-dev libz-dev 

# Install fastq-tools
RUN wget http://homes.cs.washington.edu/~dcjones/fastq-tools/fastq-tools-0.8.tar.gz
RUN tar xvf fastq-tools-0.8.tar.gz
WORKDIR fastq-tools-0.8
RUN ./configure
RUN  make install fastq==0.8

WORKDIR /work

# Define the iRODS package.
ENV ICMD_BASE="https://files.renci.org/pub/irods/releases/4.1.10/ubuntu14"
ENV ICMD_PKG="irods-icommands-4.1.10-ubuntu14-x86_64.deb"

# Install icommands.
RUN curl -o "$ICMD_PKG" "$ICMD_BASE/$ICMD_PKG" \
        && dpkg -i "$ICMD_PKG" \
        && rm -f "$ICMD_PKG"

# Install the wrapper script and the script to upload the output files.
ADD wrapper /usr/bin/wrapper
ADD upload-files /usr/bin/upload-files

# Make the wrapper script the default command.
CMD ["wrapper"]
