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
RUN chmod +x /usr/bin/Sailfish_align.pl
RUN chmod +x /usr/bin/Sailfish_align_with_gz.pl  
RUN git clone https://github.com/kingsfordgroup/sailfish.git

WORKDIR /sailfish
RUN mkdir build

WORKDIR build
RUN cmake -DFETCH_BOOST=TRUE ..
RUN make
RUN make install 
RUN export PATH=/sailfish/bin/:$PATH
RUN export LD_LIBRARY_PATH=/sailfish/lib/:$PATH #\
#&& set +o posix

WORKDIR /
ENTRYPOINT ["sailfish/bin/sailfish"]
CMD ["-h"]
