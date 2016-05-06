FROM ubuntu:14.04

MAINTAINER Chris Friedline <cfriedline@vcu.edu> # Modified by Upendra Devisetty

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


ENV BINPATH /usr/bin

RUN \
    wget http://busco.ezlab.org/files/BUSCO_v1.1b1.tar.gz && \
    tar zxvf BUSCO_v1.1b1.tar.gz

WORKDIR /BUSCO_v1.1b1

RUN \
    #sed 's|/bin/python|usr/bin/python3|' BUSCO_v1.1b1.py > temp && mv temp BUSCO_v1.1b1.py \
    chmod +x BUSCO_v1.1b1.py && cp BUSCO_v1.1b1.py $BINPATH

RUN \
    wget ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.2.30/ncbi-blast-2.2.30+-x64-linux.tar.gz && \
    wget http://eddylab.org/software/hmmer3/3.1b2/hmmer-3.1b2-linux-intel-x86_64.tar.gz && \
    wget http://bioinf.uni-greifswald.de/augustus/binaries/augustus-3.2.1.tar.gz

RUN \
    wget http://busco.ezlab.org/files/arthropoda_buscos.tar.gz && \
    wget http://busco.ezlab.org/files/vertebrata_buscos.tar.gz && \
    wget http://busco.ezlab.org/files/fungi_buscos.tar.gz && \
    wget http://busco.ezlab.org/files/bacteria_buscos.tar.gz && \
    wget http://busco.ezlab.org/files/metazoa_buscos.tar.gz && \
    wget http://busco.ezlab.org/files/eukaryota_buscos.tar.gz

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


ENV AUGUSTUS_CONFIG_PATH=/BUSCO_v1.1b1/augustus-3.2.1/config/
ENV PATH=/BUSCO_v1.1b1/ncbi-blast-2.2.30+/bin/:$PATH
ENV PATH=/BUSCO_v1.1b1/augustus-3.2.1/bin/:$PATH

ENTRYPOINT ["python3", "/usr/bin/BUSCO_v1.1b1.py"]
CMD ["-h"]
