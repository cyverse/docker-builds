## Dockerfile
FROM ubuntu:18.04
MAINTAINER Amanda Cooksey	
LABEL Description="snphylo"

ENV DEBIAN_FRONTEND=noninteractive

# Install all the updates and download dependencies
RUN apt-get update && \
    apt-get install -y \
    git \
    wget \
    bzip2 \
    python \
    nano \
    r-base \
    r-base-dev \
    r-cran-rgl \
    curl \
    unzip \
    libssl-dev \
    libcurl4-openssl-dev \
    libxml2-dev



RUN R -e "install.packages(c('getopt', 'BiocManager'), dependencies=TRUE,  repos='http://cran.rstudio.com/')"

RUN R -e "BiocManager::install(c('SNPRelate', 'gdsfmt', 'phangorn'))"

WORKDIR /usr/local/bin/

RUN git clone https://github.com/thlee/SNPhylo.git

ENV SNPHYLO_HOME="/usr/local/bin/SNPhylo/"

RUN mkdir -p SNPhylo/bin


RUN curl -O http://www.drive5.com/muscle/downloads3.8.31/muscle3.8.31_i86linux64.tar.gz \
&& tar xf muscle3.8.31_i86linux64.tar.gz -C /usr/local/bin/SNPhylo/bin \
&& ln -sf /usr/local/bin/SNPhylo/bin/muscle3.8.31_i86linux64 /usr/local/bin/SNPhylo/bin/muscle

WORKDIR /usr/local/bin/SNPhylo/bin/

RUN curl -O https://evolution.gs.washington.edu/phylip/download/phylip-3.695.zip \
&& unzip phylip-3.695.zip \
&& ln -sf /usr/local/bin/SNPhylo/bin/phylip-3.695 /usr/local/bin/SNPhylo/bin/phylip \
&& cd /usr/local/bin/SNPhylo/bin/phylip/src \
&& rm -r /usr/local/bin/SNPhylo/bin/phylip/exe/testdata \
&& cp Makefile.unx Makefile \
&& make install


ENV PATH=$PATH:/usr/local/bin/SNPhylo/bin/:/usr/local/bin/SNPhylo:/usr/local/bin/SNPhylo/bin/phylip/exe/:/usr/local/lib/R/site-library

WORKDIR /usr/local/bin/SNPhylo/

RUN bash setup.sh

WORKDIR /workdir




