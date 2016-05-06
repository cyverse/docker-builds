FROM ubuntu:14.04

MAINTAINER Upendra Devisetty <upendra@cyverse.org> 

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
    wget http://buscos.ezlab.org/files/plant_early_release.tar.gz && \
    tar xvf plant_early_release.tar.gz

WORKDIR /plant_early_release

RUN \
    chmod +x BUSCO_plants.py && cp BUSCO_plants.py $BINPATH

RUN \
    wget ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.2.30/ncbi-blast-2.2.30+-x64-linux.tar.gz && \
    wget http://eddylab.org/software/hmmer3/3.1b2/hmmer-3.1b2-linux-intel-x86_64.tar.gz && \
    wget http://bioinf.uni-greifswald.de/augustus/binaries/augustus-3.2.1.tar.gz

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


ENV AUGUSTUS_CONFIG_PATH=/plant_early_release/augustus-3.2.1/config/
ENV PATH=/plant_early_release/ncbi-blast-2.2.30+/bin/:$PATH
ENV PATH=/plant_early_release/augustus-3.2.1/bin/:$PATH

ENTRYPOINT ["python3", "/usr/bin/BUSCO_plants.py"]
CMD ["-h"]

# Building image
# docker build -t"=ubuntu/busco:2.0" .
# Runing the executables of the image without any arguments
# docker run ubuntu/busco:2.0
# Running the image with arguments
# First example dataset
# sudo docker run --rm -v $(pwd):/home/upendra_35/BUSCO/plant_early_release/sample_data -w /home/upendra_35/BUSCO/plant_early_release/sample_data ubuntu/busco:2.0 -in target.fa -l example -m genome -o SAMPLE_example -f
# Plantae dataset
# sudo docker run --rm -v $(pwd):/home/upendra_35/BUSCO/plant_early_release/sample_data -w /home/upendra_35/BUSCO/plant_early_release/sample_data ubuntu/busco:2.0 -in target.fa -l plantae -m genome -o SAMPLE_plantae -f

