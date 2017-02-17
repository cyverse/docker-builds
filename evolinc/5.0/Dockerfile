FROM ubuntu:14.04.3
MAINTAINER Upendra Devisetty <upendra@cyverse.org>
LABEL Description "This Dockerfile is for evolinc-ii pipeline"

RUN apt-get update && apt-get install -y g++ \
		make \
		git \
		zlib1g-dev \
		python \
		perl \
		wget \
		curl \
		python-matplotlib \
		python-numpy \
                python-pandas

ENV BINPATH /usr/bin
WORKDIR /evolinc_docker

# Cufflinks
RUN wget -O- http://cole-trapnell-lab.github.io/cufflinks/assets/downloads/cufflinks-2.2.1.Linux_x86_64.tar.gz | tar xzvf -

# Transdecoder
RUN wget -O- https://github.com/TransDecoder/TransDecoder/archive/2.0.1.tar.gz | tar xzvf -

# NCBI Blast
RUN curl ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.4.0/ncbi-blast-2.4.0+-x64-linux.tar.gz > ncbi-blast-2.4.0+-x64-linux.tar.gz
RUN tar xvf ncbi-blast-2.4.0+-x64-linux.tar.gz

# Samtools
RUN wget --no-check-certificate http://sourceforge.net/projects/samtools/files/samtools/1.0/samtools-bcftools-htslib-1.0_x64-linux.tar.bz2/download
RUN tar xvf download

# Bedtools
RUN wget https://github.com/arq5x/bedtools2/archive/v2.25.0.tar.gz
RUN tar xvf v2.25.0.tar.gz
RUN cd bedtools2-2.25.0 && make
RUN cd ..

# Bedops tool
RUN wget -O- https://github.com/bedops/bedops/releases/download/v2.4.16/bedops_linux_x86_64-v2.4.16.tar.bz2 | tar jxvf -

# cpan
RUN curl -L http://cpanmin.us | perl - App::cpanminus
RUN cpanm URI/Escape.pm

# R libraries
RUN echo "deb http://cran.cnr.berkeley.edu/bin/linux/ubuntu trusty/" >> /etc/apt/sources.list
RUN apt-key adv --keyserver keyserver.ubuntu.com --recv-keys 51716619E084DAB9
RUN apt-get update
RUN apt-get install -y r-base r-base-dev
RUN Rscript -e 'install.packages("splitstackshape", dependencies = TRUE, repos="http://cran.rstudio.com/");'
RUN Rscript -e 'install.packages("dplyr", dependencies = TRUE, repos="http://cran.rstudio.com/");'
RUN Rscript -e 'source("https://bioconductor.org/biocLite.R"); biocLite("Biostrings");'
RUN Rscript -e 'install.packages("getopt", dependencies = TRUE, repos="http://cran.rstudio.com/");'

# Uniprot database
ADD https://github.com/iPlantCollaborativeOpenSource/docker-builds/releases/download/evolinc-I/uniprot_sprot.tar.gz /evolinc_docker/
RUN tar zxvf /evolinc_docker/uniprot_sprot.tar.gz && rm /evolinc_docker/uniprot_sprot.tar.gz

# Biopython
RUN curl "https://bootstrap.pypa.io/get-pip.py" -o "get-pip.py"
RUN python get-pip.py
RUN pip install biopython

# Evolinc wrapper scripts
ADD *.sh *.py *.R /evolinc_docker/
RUN chmod +x /evolinc_docker/evolinc-part-I.sh && cp /evolinc_docker/evolinc-part-I.sh $BINPATH

WORKDIR /

# Setting paths to all the softwares
ENV PATH /evolinc_docker/cufflinks-2.2.1.Linux_x86_64/:$PATH
ENV PATH /evolinc_docker/TransDecoder-2.0.1/:$PATH
ENV PATH /evolinc_docker/ncbi-blast-2.4.0+/bin/:$PATH
ENV PATH /evolinc_docker/bedtools2-2.25.0/bin/:$PATH
ENV PATH /evolinc_docker/samtools-bcftools-htslib-1.0_x64-linux/bin/:$PATH
ENV PATH /evolinc_docker/bin/:$PATH

# Entrypoint
ENTRYPOINT ["evolinc-part-I.sh"]
CMD ["-h"]
