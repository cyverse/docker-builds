# Dockerfile for RseqFilt container
# load ubuntu 18.04
FROM ubuntu:18.04

MAINTAINER Renesh Bedre (reneshbe@gmail.com)

# update and install essential modules
# container timezone
ENV TZ=America/Chicago
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone
RUN apt-get update && apt-get install -y \
    git \
    python3 \
    python3-pip \
    python3-tk \
    tzdata
RUN pip3 install --upgrade pip setuptools
RUN pip3 install \
    numpy \
    termcolor \
    pysam \
    datetime \
    matplotlib

# set ENV variables and working directory
ENV PROJDIR /project
RUN mkdir -p $PROJDIR
WORKDIR $PROJDIR

# clone and install RseqFilt
RUN git clone https://github.com/reneshbedre/RseqFilt.git
WORKDIR $PROJDIR/RseqFilt
RUN python3 setup.py install

# execute filter.py command when container runs
ENTRYPOINT ["filter.py"]
