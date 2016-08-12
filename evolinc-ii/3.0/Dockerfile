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

# NCBI Blast download
RUN curl ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.3.0/ncbi-blast-2.3.0+-x64-linux.tar.gz > ncbi-blast-2.3.0+-x64-linux.tar.gz
RUN tar xvf ncbi-blast-2.3.0+-x64-linux.tar.gz

# Bedtools
RUN wget https://github.com/arq5x/bedtools2/releases/download/v2.26.0/bedtools-2.26.0.tar.gz
RUN tar zxvf bedtools-2.26.0.tar.gz
RUN cd bedtools2 && make
RUN cd ..

# Cufflinks
RUN wget -O- http://cole-trapnell-lab.github.io/cufflinks/assets/downloads/cufflinks-2.2.1.Linux_x86_64.tar.gz | tar xzvf -

# Mafft
RUN apt-get install -y mafft

# cpan
RUN apt-get install -y cpanminus

# Install BioPerl dependancies, mostly from cpan
RUN apt-get install --yes \
 libpixman-1-0 \
 libpixman-1-dev \
 graphviz \
 libxml-parser-perl \
 libsoap-lite-perl 

RUN cpanm Test::Most \
 Algorithm::Munkres \
 Array::Compare Clone \
 PostScript::TextBlock \
 SVG \
 SVG::Graph \
 Set::Scalar \
 Sort::Naturally \
 Graph \
 GraphViz \
 HTML::TableExtract \
 Convert::Binary::C \
 Math::Random \
 Error \
 Spreadsheet::ParseExcel \
 XML::Parser::PerlSAX \
 XML::SAX::Writer \
 XML::Twig XML::Writer

RUN apt-get install -y \
 libxml-libxml-perl \
 libxml-dom-xpath-perl \
 libxml-libxml-simple-perl \
 libxml-dom-perl

# Install BioPerl last built
RUN cpanm -v  \
 CJFIELDS/BioPerl-1.6.924.tar.gz 

# Setting paths to all the softwares
ENV PATH /ncbi-blast-2.3.0+/bin/:$PATH
ENV PATH /bedtools2/bin/:$PATH
ENV PATH /cufflinks-2.2.1.Linux_x86_64/:$PATH

# Biopython
RUN curl "https://bootstrap.pypa.io/get-pip.py" -o "get-pip.py"
RUN python get-pip.py
RUN pip install biopython

# R libraries
RUN echo "deb http://cran.cnr.berkeley.edu/bin/linux/ubuntu trusty/" >> /etc/apt/sources.list
RUN apt-key adv --keyserver keyserver.ubuntu.com --recv-keys 51716619E084DAB9
RUN apt-get update
RUN apt-get install -y r-base r-base-dev
RUN Rscript -e 'install.packages("getopt", dependencies = TRUE, repos="http://cran.rstudio.com/");'
RUN Rscript -e 'install.packages("reshape2", dependencies = TRUE, repos="http://cran.rstudio.com/");'

# Add all the scripts to the root directory PATHSSSS
ADD *.py *.pl *.R *.sh /
RUN chmod +x /Building_Families.sh
RUN chmod +x /evolinc-part-II.sh && cp /evolinc-part-II.sh $BINPATH

ENTRYPOINT ["/evolinc-part-II.sh"]
CMD ["-h"]

# Docker build
# docker build -t"=ubuntu/evolinc-ii:0.1" .
# Run without any commands
# sudo docker run ubuntu/evolinc-ii:0.1 -h
# Testing
# sudo docker run --rm -v $(pwd):/working-dir -w /working-dir ubuntu/evolinc-ii:0.1 -b sample.data/Blasting_list.txt -l sample.data/Hsap_lincRNAs.fasta -q Hsap -i sample.data -s sample.data/Species_list.txt -o test_out
# docker build -t"=ubuntu/evolinc-ii:0.3" .
#docker run --rm -v $(pwd):/working-dir -w /working-dir ubuntu/evolinc-ii:0.3 -b sample.data/Blasting_list.txt -l sample.data/Hsap_lincRNAs.fasta -q Hsap -i sample.data -s sample.data/Species_list.txt -e 1e-20 -o test_out

