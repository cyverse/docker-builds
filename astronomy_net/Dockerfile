#Pull base image
FROM python:2.7

MAINTAINER Andre Mercer <amercer@iplantcollaborative.org>

# Run updates
RUN apt-get update
RUN apt-get upgrade -y
RUN apt-get install wget -y

# Install dependencies
RUN apt-get upgrade gcc
RUN apt-get install -y libcairo2-dev libnetpbm10-dev netpbm \
                       libpng12-dev libjpeg-dev python-numpy \
                       python-pyfits python-dev zlib1g-dev \
                       libbz2-dev swig libcfitsio3-dev \
                       python-setuptools python-scipy \ 
                       libatlas-base-dev gfortran build-essential

RUN pip install scipy
RUN pip install numpy

RUN mkdir /usr/include/netpbm
RUN ln -s /usr/include/pam.h /usr/include/netpbm/pam.h
RUN ln -s /usr/include/ppm.h /usr/include/netpbm/ppm.h

RUN wget http://astrometry.net/downloads/astrometry.net-latest.tar.gz 
RUN tar zxf astrometry.net-latest.tar.gz
RUN make -C astrometry.net-0.56
RUN make py -C astrometry.net-0.56
RUN make extra -C astrometry.net-0.56
RUN make install -C astrometry.net-0.56

env PATH /usr/local/astrometry/bin/:$PATH
