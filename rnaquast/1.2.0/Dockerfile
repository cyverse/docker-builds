FROM ubuntu:14.04.3
MAINTAINER Upendra Devisetty <upendra@cyverse.org>
LABEL Description "This Dockerfile is for rnaQUAST-1.2.0"

# Install the dependencies
RUN apt-get update 
RUN apt-get install -y wget python-pip python-matplotlib make unzip samtools emboss emboss-lib zlib1g-dev
RUN pip install gffutils 
RUN pip install joblib 

# Install rnaQUAST-1.2.0
RUN wget -O- http://spades.bioinf.spbau.ru/rnaquast/release1.2.0/rnaQUAST-1.2.0.tar.gz | tar zxvf -
RUN chmod +x /rnaQUAST-1.2.0/rnaQUAST.py

# Set working directory
WORKDIR /rnaQUAST-1.2.0

# Install other softwares
# Blast
RUN wget -O- http://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ncbi-blast-2.3.0+-x64-linux.tar.gz | tar zxvf  -

# GMAP (The executable is already in path)
RUN wget -O- http://research-pub.gene.com/gmap/src/gmap-gsnap-2015-12-31.v6.tar.gz | tar zxvf -
RUN cd gmap-2015-12-31 && ./configure && make && make check && make install

# BLAT
RUN mkdir Blat
WORKDIR Blat
RUN wget http://hgwdev.cse.ucsc.edu/~kent/exe/linux/blatSuite.36.zip
RUN unzip blatSuite.36.zip
RUN chmod +x blat

# TopHat read alignment
WORKDIR /rnaQUAST-1.2.0
RUN wget -O- https://ccb.jhu.edu/software/tophat/downloads/tophat-2.1.0.Linux_x86_64.tar.gz | tar zxvf -

# Bowtie2
WORKDIR /rnaQUAST-1.2.0
RUN wget https://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.2.7/bowtie2-2.2.7-linux-x86_64.zip/download
RUN unzip download

# BUSCO
RUN mkdir BUSCO_v1.1b1 
WORKDIR BUSCO_v1.1b1
RUN wget https://raw.githubusercontent.com/upendrak/rnaQUAST-1.2.0/master/BUSCO_v1.1b1.py
RUN chmod +x BUSCO_v1.1b1.py

# Hmmer: (it should be set to path)
WORKDIR /rnaQUAST-1.2.0
RUN wget -O- http://eddylab.org/software/hmmer3/3.1b2/hmmer-3.1b2-linux-intel-x86_64.tar.gz | tar zxvf -
RUN cd hmmer-3.1b2-linux-intel-x86_64 && ./configure && make && make check && make install    

# GeneMarkS-T
RUN mkdir GeneMarkS-T
WORKDIR GeneMarkS-T
RUN wget -O- http://topaz.gatech.edu/GeneMark/tmp/GMtool_BJnvL/gmst_linux_64.tar.gz | tar zxvf -
RUN chmod +x gmst.pl

# STAR for read alignment
WORKDIR /rnaQUAST-1.2.0
RUN wget -O- https://github.com/alexdobin/STAR/archive/2.5.1b.tar.gz | tar zxvf -
RUN cd STAR-2.5.1b && make STAR

# Set environmental paths
ENV PATH /rnaQUAST-1.2.0/Blat/:$PATH
ENV PATH /rnaQUAST-1.2.0/ncbi-blast-2.3.0+/bin/:$PATH
ENV PATH /rnaQUAST-1.2.0/tophat-2.1.0.Linux_x86_64/:$PATH
ENV PATH /rnaQUAST-1.2.0/bowtie2-2.2.7/:$PATH
ENV PATH /rnaQUAST-1.2.0/BUSCO_v1.1b1/:$PATH
ENV PATH /rnaQUAST-1.2.0/GeneMarkS-T/:$PATH
ENV PATH /rnaQUAST-1.2.0/STAR-2.5.1b/bin/Linux_x86_64_static/:$PATH

ENTRYPOINT ["/rnaQUAST-1.2.0/rnaQUAST.py"]
CMD ["-h"]
