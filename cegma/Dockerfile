## Dockerfile
FROM ubuntu:16.04
MAINTAINER Amanda Cooksey	
LABEL Description="CEGMA no longer supported. use BUSCO"

# Install all the updates and download dependencies
RUN apt-get update
RUN apt-get install --yes \
 wget \
 gcc \
 nano \
 apt-utils \
 build-essential \
 perl 

RUN cpan inc::latest Getopt::Std

RUN echo 'export PATH=/opt/conda/bin:$PATH' > /etc/profile.d/conda.sh && \
    wget --quiet https://repo.continuum.io/miniconda/Miniconda2-4.0.5-Linux-x86_64.sh -O ~/miniconda.sh && \
    /bin/bash ~/miniconda.sh -b -p /opt/conda && \
    rm ~/miniconda.sh

RUN chmod 777 -R /opt/conda/

ENV PATH=$PATH:/opt/conda/bin

RUN conda config --add channels bioconda

RUN conda upgrade conda

RUN conda install -c conda-forge -c bioconda wise2==2.4.1

ADD geneid_v1.4.4.Jan_13_2011.tar.gz ncbi-blast-2.2.27+-x64-linux.tar.gz hmmer-3.0-linux-intel-x86_64.tar.gz /usr/bin/

ADD CEGMA_v2.5.tar.gz /work-dir

RUN cd /work-dir/CEGMA_v2.5 && make

RUN cd /usr/bin/geneid && make 

RUN cd /usr/bin/hmmer-3.0-linux-intel-x86_64 && ./configure && make && make install

RUN cd /work-dir/CEGMA_v2.5 && make

ENV CEGMA="/work-dir/CEGMA_v2.5"
ENV CEGMATMP="/work-dir"
ENV PERL5LIB="$PERL5LIB:$CEGMA/lib"
ENV WISECONFIGDIR="/opt/conda/share/wise2/wisecfg"

ENV PATH=$PATH:/usr/bin/ncbi-blast-2.2.27+/bin:/usr/bin/hmmer-3.0-linux-intel-x86_64/binaries:/usr/bin/geneid/bin:/work-dir/CEGMA_v2.5/bin

# Entrypoint
ENTRYPOINT ["cegma"]

WORKDIR /work-dir
