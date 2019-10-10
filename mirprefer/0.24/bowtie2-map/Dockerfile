FROM upendradevisetty/bowtie2:2.3.4.1
MAINTAINER Upendra Kumar Devisetty <upendra@cyverse.org>
LABEL Description "This Dockerfile is for bowite2-miRPREFeR"

ENV LANG=C.UTF-8 LC_ALL=C.UTF-8

RUN apt-get update --fix-missing && \
    apt-get install -y wget \
		       unzip \
    		       liblzma-dev \
                       libbz2-dev \
                       zlib1g-dev \
                       g++ \
		       bzip2 \
		       make \
		       python 

WORKDIR /
RUN wget https://github.com/samtools/samtools/releases/download/1.8/samtools-1.8.tar.bz2 && \
    bzip2 -d samtools-1.8.tar.bz2 && \
    tar xvf samtools-1.8.tar && \
    cd samtools-1.8 && \
    ./configure --without-curses && \
    make && \
    make install

ADD bowtie2_wrapper.sh /usr/bin/
RUN chmod +x /usr/bin/bowtie2_wrapper.sh

ENTRYPOINT ["bowtie2_wrapper.sh"]
CMD [ "-h" ]
