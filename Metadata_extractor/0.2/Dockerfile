From ubuntu:16.04
MAINTAINER HUANIAN ZHANG <fantasyzhn@email.arizona.edu>
LABEL "This Dockerfile is for Metadata_extractor:ME_v02"

RUN apt-get update && apt-get install -y git make gcc wget python-dev python2.7 python-numpy python-pip
RUN pip install --no-deps astropy 
RUN pip install fitsio==0.9.7
RUN pip install matplotlib
RUN apt-get install -y python-tk
ENV BINPATH /usr/bin
ADD extract_metadata.py $BINPATH
RUN chmod +x $BINPATH/extract_metadata.py

ENTRYPOINT ["extract_metadata.py"]
