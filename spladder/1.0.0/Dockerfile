FROM ubuntu:14.04.3
MAINTAINER Upendra Devisetty <upendra@cyverse.org>
LABEL Description "This Dockerfile is for Spladder 1.0.0"

RUN apt-get update
RUN apt-get install --yes build-essential python2.7-dev python-numpy python-matplotlib python-pip zlib1g-dev libhdf5-dev git libblas-dev liblapack-dev libatlas-base-dev gfortran wget
RUN pip install pysam==0.9.0 h5py Cython git+https://github.com/scipy/scipy.git
# Spladder
RUN wget https://github.com/ratschlab/spladder/archive/v1.0.0.tar.gz
RUN tar xvf v1.0.0.tar.gz
RUN cp spladder-1.0.0/python/spladder.py /usr/bin
RUN chmod +x /usr/bin/spladder.py
# Wrapper script for Spladder
ADD spladder_wrapper_script.sh /usr/bin
RUN chmod +x /usr/bin/spladder_wrapper_script.sh
ENTRYPOINT ["/usr/bin/spladder_wrapper_script.sh"]
