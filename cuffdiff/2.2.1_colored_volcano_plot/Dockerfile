FROM ubuntu:16.04
MAINTAINER Upendra Devisetty <upendra@cyverse.org>
LABEL version="2.2.1" description='This Dockerfile is a wrapper for Cuffdiff2-16-way-max app'

# Install the updates and dependencies
RUN apt-get update && apt-get install -y wget g++ make r-base r-bioc-cummerbund

# Dumper
RUN wget http://www.cpan.org/modules/by-module/Data/Data-Dumper-2.161.tar.gz && \
	tar xvf Data-Dumper-2.161.tar.gz
WORKDIR Data-Dumper-2.161
RUN perl Makefile.PL && make && make install

# Download the cufflinks software
WORKDIR /
COPY cufflinks-1.3.0.Linux_x86_64.tar.gz /
RUN tar xvf /cufflinks-1.3.0.Linux_x86_64.tar.gz
RUN wget http://cole-trapnell-lab.github.io/cufflinks/assets/downloads/cufflinks-2.0.2.Linux_x86_64.tar.gz && \
	tar xvf /cufflinks-2.0.2.Linux_x86_64.tar.gz
RUN wget http://cole-trapnell-lab.github.io/cufflinks/assets/downloads/cufflinks-2.1.1.Linux_x86_64.tar.gz && \
	tar xvf /cufflinks-2.1.1.Linux_x86_64.tar.gz
RUN wget http://cole-trapnell-lab.github.io/cufflinks/assets/downloads/cufflinks-2.2.0.Linux_x86_64.tar.gz && \
	tar xvf /cufflinks-2.2.0.Linux_x86_64.tar.gz
RUN wget http://cole-trapnell-lab.github.io/cufflinks/assets/downloads/cufflinks-2.2.1.Linux_x86_64.tar.gz && \
	tar xvf /cufflinks-2.2.1.Linux_x86_64.tar.gz

# Copy the scripts and folders
COPY cuffdiff_fix_annotations.pl cuffdiff_get_goodstuff_test.pl cuffdiff_test16.pl munge_ids.pl /
ADD ensembl_plant_gene_desc.tar.gz /

# ENTRYPOINT
ENTRYPOINT ["/cuffdiff_test16.pl"]

