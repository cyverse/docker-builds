FROM ubuntu:14.04.3
MAINTAINER Upendra Devisetty <upendra@cyverse.org>
LABEL Description "This Dockerfile is for tbl2asn for converting fasta to .sqn format"

RUN apt-get update 

# Dependencies
RUN apt-get install -y wget
                          
# Download and change permissions of tbl2asn software
RUN wget ftp://ftp.ncbi.nih.gov/toolbox/ncbi_tools/converters/by_program/tbl2asn/linux64.tbl2asn.gz
RUN gunzip linux64.tbl2asn.gz
RUN mv linux64.tbl2asn tbl2asn
RUN chmod a+x tbl2asn
RUN cp tbl2asn /usr/bin

# Specify entrypoint
ENTRYPOINT ["tbl2asn"]
CMD ["-h"]

# Build the image
# docker build -t"=ubuntu/tbl2asn" .
# Testing the image
# Without any arguments
# sudo docker run docker run ubuntu/tbl2asn -h
# With arguments
# All the files in the current directory
# sudo docker run --rm -v $(pwd):/working-dir -w /working-dir ubuntu/tbl2asn -t template.sbt -p . -a s -V v -j "[organism=Saccharomyces cerevisiae] [strain=S288C]" 

# Testing in a directory of files
# mkdir dir && cp test.fsa dir
# sudo docker run --rm -v $(pwd):/working-dir -w /working-dir ubuntu/tbl2asn -t template.sbt -p dir -a s -V v -j "[organism=Saccharomyces cerevisiae] [strain=S288C]"

