# Sailfish 
# VERSION               0.9.2
# 
FROM      ubuntu:14.04.3 
MAINTAINER Kapeel Chougule
LABEL Description="This image is used for running Sailfish RNA seq qauntification tool "
RUN apt-get update && apt-get install -y build-essential git cmake curl unzip zlib1g-dev wget python libtbb-dev libjemalloc-dev libgoogle-perftools-dev \
pkg-config libbz2-dev 
#libhdf5-dev pkg-config liblzma-dev apt-get install libbz2-dev curl unzip
ADD Sailfish_align.pl /usr/bin/ 
ADD Sailfish_align_with_gz.pl /usr/bin/
RUN [ "chmod", "+x",  "/usr/bin/Sailfish_align.pl" ]
RUN [ "chmod", "+x",  "/usr/bin/Sailfish_align_with_gz.pl" ] \  
&& git clone https://github.com/kingsfordgroup/sailfish.git \
&& cd sailfish \
&& mkdir build \
&& cd build \
&& cmake -DFETCH_BOOST=TRUE .. \
&& make \
&& make install \
&& export PATH=$PATH:/root/sailfish/bin \
&& set +o posix
#ENTRYPOINT ["sailfish"]
#CMD ["-h"]
