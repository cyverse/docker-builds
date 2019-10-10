################################################################################
##################### Set Inital Image to work from ############################
FROM ubuntu:18.04

################################################################################
##################### Add Container Labels #####################################
LABEL "Regtools_License"="MIT"
LABEL "Description"="Software package which integrate DNA-seq and RNA-seq data\
                     to help interpret mutations in a regulatory and splicing\
                     context."
MAINTAINER Upendra Devisetty <upendra@cyverse.org>

################################################################################
##################### Install System Dependencies ##############################
RUN apt-get update -y && apt-get install -y \
    git \
    build-essential \
    cmake \
    zlib1g-dev
RUN apt-get install -y wget
################################################################################
##################### Install Regtools #########################################

# clone git repository
RUN wget https://github.com/griffithlab/regtools/archive/0.5.0.tar.gz
RUN tar xvf 0.5.0.tar.gz

# make a build directory for regtools
WORKDIR /regtools-0.5.0
RUN mkdir build
WORKDIR /regtools-0.5.0/build

# compile from source
RUN cmake .. && make
################################################################################
###################### set environment path    #################################

# add regtools executable to path
ENV PATH="/regtools-0.5.0/build:${PATH}"


# Entrypoint
ENTRYPOINT ["regtools"]
