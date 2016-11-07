FROM ubuntu:14.04.3
MAINTAINER Andre Mercer <amercer1@email.arizona.edu>
LABEL Description "This Dockerfile is for spladder 1.0"

# Install the dependencies
RUN apt-get update 
RUN apt-get install -y wget python-pip python-dev python-setuptools python-matplotlib make unzip samtools emboss emboss-lib zlib1g-dev
RUN apt-get install -y python-virtualenv git python-software-properties software-properties-common libcurl4-openssl-dev
RUN apt-get install -y libpng-dev libfreetype6-dev
RUN apt-get install -y python-scipy 
RUN apt-get install -y python-numpy

## setup virtualenv
RUN virtualenv VIRTUAL
RUN VIRTUAL/bin/pip install --upgrade pip
RUN VIRTUAL/bin/pip install scipy 
RUN VIRTUAL/bin/pip install numpy
RUN VIRTUAL/bin/pip install pysam
RUN VIRTUAL/bin/pip install h5py
RUN VIRTUAL/bin/pip install matplotlib

## Install spladder
RUN git clone https://github.com/ratschlab/spladder.git
ENTRYPOINT ["/VIRTUAL/bin/python", "/spladder/python/spladder.py"]
