FROM ubuntu:14.04.3
MAINTAINER Benjamin Bolduc <bolduc.10@osu.edu>

# Remove those annoying warnings
ENV DEBIAN_FRONTEND noninteractive

ENV BINPATH /usr/bin

RUN apt-get update && apt-get install -y \
    automake \
    build-essential \
    zlib1g-dev \
    libncurses5-dev \
    wget \
    unzip \
    git \
    python \
    python-numpy \
    python-dev \
    python-pandas \
    python-matplotlib \
    libtool
    
# Install bwa
RUN wget --no-verbose https://github.com/lh3/bwa/releases/download/v0.7.13/bwa-0.7.13.tar.bz2
RUN tar xf bwa-0.7.13.tar.bz2 && cd bwa-0.7.13 && make && cp bwa $BINPATH

# Install Bowtie2
RUN wget --no-verbose http://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.2.6/bowtie2-2.2.6-linux-x86_64.zip
RUN unzip bowtie2-2.2.6-linux-x86_64.zip && cp bowtie2-2.2.6/bowtie2* $BINPATH

# Install Samtools
RUN wget --no-verbose https://github.com/samtools/samtools/releases/download/1.3/samtools-1.3.tar.bz2
RUN tar xf samtools-1.3.tar.bz2 && cd samtools-1.3 && make && make install

# What you get for installing specific version of python through apt-get
#RUN ln -s /usr/bin/python2.7 /usr/bin/python

# Install pip
RUN wget --no-verbose https://bootstrap.pypa.io/get-pip.py
RUN python get-pip.py

# InsecurePlatformWarning for pip running under python < 2.7.9
# https://github.com/pypa/pip/issues/2681

# Install pysam (not a dependency of BamM... but needed for running tests?)
RUN pip install pysam

# Install seaborn and palettable for generating coverage figure
RUN pip install seaborn
RUN pip install palettable

# Install BamM
RUN wget --no-verbose -O BamM-v1.7.0.tar.gz https://github.com/Ecogenomics/BamM/archive/v1.7.0.tar.gz
RUN tar xf BamM-v1.7.0.tar.gz && cd BamM-1.7.0 && python setup.py install

# Clean stuff up
RUN apt-get clean && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

COPY scripts/ /usr/bin/
RUN chmod +x /usr/bin/*.py

# # Clean up fonts so matplotlib can find them
RUN rm -rf ~/.cache/matplotlib

ENTRYPOINT ["Read2ReferenceMapper.py"]
CMD ["--help"]