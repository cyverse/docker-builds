FROM ubuntu:focal
LABEL Description "This Dockerfile is for blast-2.13.0"

RUN apt-get update && apt-get install  -y \
		wget \
		python3 \
		libgomp1 \
	&& rm -rf /var/lib/apt/lists/*

# Download & install BLAST
RUN mkdir /opt/blast
WORKDIR /opt/blast
RUN wget https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ncbi-blast-2.13.0+-x64-linux.tar.gz
RUN tar -zxvf /opt/blast/ncbi-blast-2.13.0+-x64-linux.tar.gz
ENV PATH /opt/blast/ncbi-blast-2.13.0+/bin/:$PATH
WORKDIR /

# Add wrapper script
ADD blastn_wrapper.sh blastp_wrapper.sh /usr/bin/
RUN chmod +x /usr/bin/blastn_wrapper.sh /usr/bin/blastp_wrapper.sh

