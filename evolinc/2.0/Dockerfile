FROM ubuntu:14.04.3
MAINTAINER Upendra Devisetty

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

# Quast
RUN wget -O- https://downloads.sourceforge.net/project/quast/quast-4.0.tar.gz | tar zxvf -

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

# Evolinc wrapper scripts
ADD *.sh *.py *.pl /evolinc_docker/
RUN chmod +x /evolinc_docker/evolinc-part-I.sh && cp /evolinc_docker/evolinc-part-I.sh $BINPATH

# Setting paths to all the softwares
ENV PATH /evolinc_docker/cufflinks-2.2.1.Linux_x86_64/:$PATH
ENV PATH /evolinc_docker/TransDecoder-2.0.1/:$PATH
ENV PATH /evolinc_docker/ncbi-blast-2.4.0+/bin/:$PATH
ENV PATH /evolinc_docker/bedtools2-2.25.0/bin/:$PATH
ENV PATH /evolinc_docker/samtools-bcftools-htslib-1.0_x64-linux/bin/:$PATH
ENV PATH /evolinc_docker/bin/:$PATH
ENV PATH /evolinc_docker/quast-4.0/:$PATH

# Entrypoint
ENTRYPOINT ["evolinc-part-I.sh"]
CMD ["-h"]

# Docker build
# docker build -t"=ubuntu/evolinc:0.2" dockerfile-evolinc
# Run it to test
# docker run -it ubuntu/evolinc:0.2
# mkdir /workind-dir
# sudo git clone https://upendra_35@bitbucket.org/upendra_35/evolinc_docker.git
# docker run --rm -v $(pwd):/working-dir -w /working-dir ubuntu/evolinc:0.2 -c AthalianaslutteandluiN30merged.gtf -g TAIR10_chr.fasta -r TAIR10_GFF3_genes_mod.gff -b TE_RNA_transcripts.fa -o test_out_new -t AnnotatedPEATPeaks.gff -x Atha_known_lncRNAs.mod.gff 
# docker tag ubuntu/evolinc:0.2 upendradevisetty/evolinc:0.2
