FROM ubuntu:16.04
MAINTAINER Upendra Devisetty <upendra@cyverse.org>
LABEL version="0.24" description="This Dockerfile is for building miRPREFER-0.24"

RUN apt-get update && apt-get install -y wget \
		      g++ \
                      bzip2 \
                      zlib1g-dev \
                      libbz2-dev \
		      make \
		      python

# Vienna Package
ENV VIENNA_VERSION=2.4.6
RUN wget https://www.tbi.univie.ac.at/RNA/download/sourcecode/2_4_x/ViennaRNA-${VIENNA_VERSION}.tar.gz
RUN tar xvf ViennaRNA-${VIENNA_VERSION}.tar.gz
WORKDIR ViennaRNA-${VIENNA_VERSION}
RUN ./configure --prefix="/ViennaRNA/" --without-perl && make && make install
WORKDIR /

# Samtools
ENV SAMTOOLS_VERSION=1.8
RUN wget https://sourceforge.net/projects/samtools/files/samtools/1.8/samtools-${SAMTOOLS_VERSION}.tar.bz2
RUN bzip2 -d samtools-${SAMTOOLS_VERSION}.tar.bz2
RUN tar xvf samtools-${SAMTOOLS_VERSION}.tar
WORKDIR samtools-${SAMTOOLS_VERSION}
RUN ./configure --without-curses --disable-lzma && make && make install
WORKDIR /

# miR-PREFER
ENV miRPREFER_VERSION=0.24
RUN wget https://github.com/hangelwen/miR-PREFeR/archive/v${miRPREFER_VERSION}.tar.gz
RUN tar xvf v${miRPREFER_VERSION}.tar.gz
COPY miR_PREFeR.py /miR-PREFeR-0.24

# Environment settings
ENV PATH /ViennaRNA/bin/:$PATH

# Wrapper script
COPY miRPREFER.sh /

# ENTRYPOINT
ENTRYPOINT ["bash", "/miRPREFER.sh"]
CMD ["-h"]
