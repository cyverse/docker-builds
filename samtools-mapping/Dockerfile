FROM ubuntu:14.04.3
MAINTAINER Upendra Devisetty <upendra@cyverse.org>
LABEL Description "This Dockerfile is for samtools_mapping_wrapper script"

RUN apt-get update && apt-get install -y samtools

ADD samtools_mapping_wrapper.sh /
RUN chmod +x /samtools_mapping_wrapper.sh && cp /samtools_mapping_wrapper.sh /usr/bin

ENTRYPOINT ["samtools_mapping_wrapper.sh"]
