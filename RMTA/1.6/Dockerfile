FROM ubuntu:14.04.3
MAINTAINER Upendra Devisetty <upendra@cyverse.org>
LABEL Description "This Dockerfile is used for hisat2 tool with sra option and cufflinks/stringtie with Cuffcompare and Cuffmerge"

RUN apt-get update && apt-get install -y build-essential \
                                         git \
                                         python \
                                         wget \
                                         unzip \
					 build-essential \
        				 zlib1g-dev \
        				 libncurses5-dev \
        				 software-properties-common

RUN add-apt-repository -y ppa:openjdk-r/ppa
RUN apt-get update
RUN apt-get install -y openjdk-8-jdk

ENV BINPATH /usr/bin

# NCBI SRA-TOOL kit
WORKDIR /ncbi
RUN git clone https://github.com/ncbi/ngs.git
RUN git clone https://github.com/ncbi/ncbi-vdb.git
RUN git clone https://github.com/ncbi/sra-tools.git
RUN ngs/ngs-sdk/configure --prefix=~/software/apps/sratoolkit/gcc/64/2.5.8
RUN make default install -C ngs/ngs-sdk
RUN ncbi-vdb/configure --prefix=~/software/apps/sratoolkit/gcc/64/2.5.8
RUN make default install -C ncbi-vdb
RUN sra-tools/configure --prefix=~/software/apps/sratoolkit/gcc/64/2.5.8
RUN make default install -C sra-tools

# Hisat2
WORKDIR /hisat2
RUN wget ftp://ftp.ccb.jhu.edu/pub/infphilo/hisat2/downloads/hisat2-2.0.5-Linux_x86_64.zip 
RUN unzip hisat2-2.0.5-Linux_x86_64.zip
RUN cp hisat2-2.0.5/hisat* $BINPATH

# Cufflinks
WORKDIR /
RUN wget -O- http://cole-trapnell-lab.github.io/cufflinks/assets/downloads/cufflinks-2.2.1.Linux_x86_64.tar.gz | tar xzvf -

# Stringtie
RUN wget -O- http://ccb.jhu.edu/software/stringtie/dl/stringtie-1.2.4.Linux_x86_64.tar.gz | tar xzvf -

# Samtools
RUN wget https://github.com/samtools/samtools/releases/download/1.3.1/samtools-1.3.1.tar.bz2 
RUN tar xvf samtools-1.3.1.tar.bz2
WORKDIR /samtools-1.3.1
RUN make

# Picard
WORKDIR /
RUN wget https://github.com/lomereiter/sambamba/releases/download/v0.6.6/sambamba_v0.6.6_linux.tar.bz2
RUN tar xvf sambamba_v0.6.6_linux.tar.bz2
RUN mv sambamba_v0.6.6 /usr/bin

# Wrapper script
ADD Hisat2-Cuffcompare-Cuffmerge.sh $BINPATH
RUN chmod +x $BINPATH/Hisat2-Cuffcompare-Cuffmerge.sh

# Set environment
ENV PATH /cufflinks-2.2.1.Linux_x86_64/:$PATH
ENV PATH /stringtie-1.2.4.Linux_x86_64/:$PATH
ENV PATH /samtools-1.3.1/:$PATH

ENTRYPOINT ["Hisat2-Cuffcompare-Cuffmerge.sh"]
