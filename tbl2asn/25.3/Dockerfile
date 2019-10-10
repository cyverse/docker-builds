FROM ubuntu:14.04.3
MAINTAINER Upendra Devisetty <upendra@cyverse.org>
LABEL Description "This Dockerfile is for tbl2asn for converting fasta to .sqn format version 25.3"

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
CMD ["--help"]
