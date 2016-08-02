# Sailfish 
# VERSION               0.9.2
# 
FROM      ubuntu:14.04.3 
MAINTAINER Upendra Kumar Devisetty <upendra@cyverse.org>
LABEL Description="This image is used for running Sailfish RNA seq qauntification tool "
RUN apt-get update && apt-get install -y build-essential git cmake curl unzip zlib1g-dev wget python libtbb-dev libjemalloc-dev libgoogle-perftools-dev \
pkg-config libbz2-dev 
#libhdf5-dev pkg-config liblzma-dev apt-get install libbz2-dev curl unzip

RUN git clone https://github.com/kingsfordgroup/sailfish.git
WORKDIR /sailfish
RUN mkdir build
WORKDIR build
RUN cmake -DFETCH_BOOST=TRUE ..
RUN make
RUN make install 
RUN export PATH=/sailfish/bin/:$PATH
RUN export LD_LIBRARY_PATH=/sailfish/lib/:$PATH 

WORKDIR /

# Add wrapper ScripT
ADD Sailfish_align.pl /usr/bin/ 
RUN chmod +x /usr/bin/Sailfish_align.pl

#ENTRYPOINT ["/usr/bin/Sailfish_align.pl"]
