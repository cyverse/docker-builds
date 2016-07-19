FROM ubuntu:14.04.3
MAINTAINER Upendra Devisetty <upendra@cyverse.org>
LABEL Description="This image is used for running BWA-0.7.15"
RUN apt-get update
RUN apt-get install --yes build-essential gcc-multilib apt-utils zlib1g-dev wget
RUN wget https://sourceforge.net/projects/bio-bwa/files/bwa-0.7.15.tar.bz2/download
RUN tar xvf download
WORKDIR /bwa-0.7.15
RUN make 
WORKDIR /
# Add wrapper script
ADD bwa_mem_nms_modified_bwa_0.7.15.pl /
RUN chmod +x /bwa_mem_nms_modified_bwa_0.7.15.pl && mv /bwa_mem_nms_modified_bwa_0.7.15.pl /usr/bin

ENTRYPOINT ["bwa_mem_nms_modified_bwa_0.7.15.pl"]
