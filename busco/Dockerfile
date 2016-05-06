FROM ubuntu:14.04

MAINTAINER Chris Friedline <cfriedline@vcu.edu>

RUN \
    sed -i 's/# \(.*multiverse$\)/\1/g' /etc/apt/sources.list && \
    sed -i 's/# \(.*universe$\)/\1/g' /etc/apt/sources.list && \
    apt-get update && \
    apt-get -y upgrade && \
    apt-get install -y curl git wget build-essential cmake

RUN \
    apt-get install -y \
    python3 \
    python3-dev \
    libboost-iostreams-dev \
    libz-dev \
    libgsl0-dev \
    libboost-graph-dev \
    samtools \
    libbam-dev \
    vim \
    emboss \
    emboss-lib

RUN \
    mkdir /busco && \
    cd /busco && \
    wget ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.2.30/ncbi-blast-2.2.30+-x64-linux.tar.gz && \
    wget http://selab.janelia.org/software/hmmer3/3.1b2/hmmer-3.1b2-linux-intel-x86_64.tar.gz && \
    wget http://bioinf.uni-greifswald.de/augustus/binaries/augustus-3.2.1.tar.gz

RUN \
    cd /busco && \
    wget http://busco.ezlab.org/files/arthropoda_buscos.tar.gz && \
    wget http://busco.ezlab.org/files/vertebrata_buscos.tar.gz && \
    wget http://busco.ezlab.org/files/fungi_buscos.tar.gz && \
    wget http://busco.ezlab.org/files/bacteria_buscos.tar.gz && \
    wget http://busco.ezlab.org/files/metazoa_buscos.tar.gz && \
    wget http://busco.ezlab.org/files/eukaryota_buscos.tar.gz

RUN \
    cd /busco && \
    for i in *.gz; do tar zxvf $i; done;

RUN \
    cd /busco/hmmer* && \
    ./configure && make && make install

RUN \
    cd /busco && \
    git clone https://github.com/pezmaster31/bamtools.git && \
    cd bamtools && mkdir build && cd build && \
    cmake .. && \
    make && \
    make install

RUN \
    ln -s /usr/local/include/bamtools /usr/include/bamtools && \
    ln -s /usr/local/lib/bamtools/* /usr/local/lib && \
    cd /busco/augustus* && \
    make

RUN \
    cd /busco && \
    wget http://busco.ezlab.org/files/BUSCO_v1.1b1.tar.gz && \
    tar zxvf BUSCO_v1.1b1.tar.gz

ENV AUGUSTUS_CONFIG_PATH=/busco/augustus-3.2.1/config/
ENV PATH=/busco/ncbi-blast-2.2.30+/bin:/busco/augustus-3.2.1/bin:$PATH

COPY BUSCO_v1.1b1.py /busco/BUSCO_v1.1b1

WORKDIR /results

CMD bash
