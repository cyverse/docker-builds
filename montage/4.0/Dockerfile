FROM ubuntu:14.04.3
MAINTAINER Huanian Zhang <fantasyzhn@email.arizona.edu>
LABEL "This Dockerfile is for Montage:4.0"

RUN apt-get update && apt-get install -y git make gcc
RUN git clone https://github.com/Caltech-IPAC/Montage
WORKDIR /Montage
#RUN cd Montage
RUN make
ENV PATH /Montage/bin:$PATH
WORKDIR /
ENTRYPOINT ["mJPEG"]

