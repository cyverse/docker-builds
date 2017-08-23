FROM ubuntu:14.04
MAINTAINER Upendra Devisetty <upendra@cyverse.org> 
LABEL Description "This Dockerfile is for BUSCO v2.0"

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
    emboss \
    emboss-lib


ENV BINPATH /usr/bin

RUN \
    wget ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.6.0/ncbi-blast-2.6.0+-x64-linux.tar.gz && \
    wget http://eddylab.org/software/hmmer3/3.1b2/hmmer-3.1b2-linux-intel-x86_64.tar.gz && \
    wget http://bioinf.uni-greifswald.de/augustus/binaries/augustus-3.3.tar.gz

RUN \
    for i in *.gz; do tar zxvf $i; done;

RUN \
    cd hmmer* && \
    ./configure && make && make install

RUN \
    git clone https://github.com/pezmaster31/bamtools.git && \
    cd bamtools && mkdir build && cd build && \
    cmake .. && \
    make && \
    make install

RUN \
    ln -s /usr/local/include/bamtools /usr/include/bamtools && \
    ln -s /usr/local/lib/bamtools/* /usr/local/lib && \
    cd augustus* && \
    make

ENV AUGUSTUS_CONFIG_PATH=/augustus/config/
ENV PATH=/augustus/bin/:$PATH
ENV PATH=/augustus/scripts/:$PATH
ENV PATH=/ncbi-blast-2.6.0+/bin/:$PATH
ENV PATH=/hmmer-3.1b2-linux-intel-x86_64/binaries/:$PATH

RUN git clone https://gitlab.com/ezlab/busco.git
WORKDIR busco
RUN python3 setup.py install
ADD config.ini config/

ENTRYPOINT ["python3", "/busco/scripts/run_BUSCO.py"]
CMD ["-h"]

