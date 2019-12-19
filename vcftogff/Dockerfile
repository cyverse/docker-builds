## Dockerfile
FROM ubuntu:16.04
MAINTAINER Amanda Cooksey	
LABEL Description="Dockerized version vcf 2 gff"

# Install all the updates and download dependencies
RUN apt-get update
RUN apt-get install --yes \
 apt-utils \
 build-essential \
 libexpat1-dev \
 perl

RUN cpan App::cpanminus

RUN cpanm --force URI::Escape Pod::Usage Getopt::Long

ADD vcf_to_gff.pl /usr/bin/

# Change the permissions and the path for the script
RUN chmod a+x /usr/bin/vcf_to_gff.pl

# Entrypoint
ENTRYPOINT ["perl", "/usr/bin/vcf_to_gff.pl"]

RUN mkdir /work-dir

WORKDIR /work-dir

