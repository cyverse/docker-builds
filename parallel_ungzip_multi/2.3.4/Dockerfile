FROM ubuntu:14.04.3
MAINTAINER Upendra Devisetty <upendra@cyverse.org>
LABEL Description "This Dockerfile is for parallel gzip"

RUN apt-get update && apt-get install -y make wget g++ zlib1g-dev

ENV BINPATH /usr/bin
ENV PIGZ_VERSION 2.3.4

# pigz software
RUN wget http://www.zlib.net/pigz/pigz-$PIGZ_VERSION.tar.gz && \
    tar xvf pigz-$PIGZ_VERSION.tar.gz && \
    rm pigz-$PIGZ_VERSION.tar.gz && \
    cd pigz-$PIGZ_VERSION && \
    make && \
    mv unpigz $BINPATH 

# Wrapper script
ADD unpigz_multi_wrapper.sh $BINPATH
RUN chmod +x $BINPATH/unpigz_multi_wrapper.sh

ENTRYPOINT ["unpigz_multi_wrapper.sh"]
CMD [-h]
