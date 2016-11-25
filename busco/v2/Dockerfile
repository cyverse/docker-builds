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
    wget https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ncbi-blast-2.5.0+-x64-linux.tar.gz && \
    wget http://eddylab.org/software/hmmer3/3.1b2/hmmer-3.1b2-linux-intel-x86_64.tar.gz && \
    wget http://bioinf.uni-greifswald.de/augustus/binaries/augustus-3.2.2.tar.gz

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

RUN wget https://gitlab.com/ezlab/busco/raw/master/BUSCO.py
RUN chmod +x \BUSCO.py && cp \BUSCO.py $BINPATH

ENV AUGUSTUS_CONFIG_PATH=/augustus-3.2.2/config/
ENV PATH=/augustus-3.2.2/bin/:$PATH
ENV PATH=/augustus-3.2.2/scripts/:$PATH
ENV PATH=/ncbi-blast-2.5.0+/bin/:$PATH

ENTRYPOINT ["python3", "/usr/bin/BUSCO.py"]
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

