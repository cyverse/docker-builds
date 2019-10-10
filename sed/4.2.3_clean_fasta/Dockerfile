FROM ubuntu:14.04.3
MAINTAINER Upendra Devisetty <upendra@cyverse.org>
LABEL Description "This Dockerfile is for sed"

RUN apt-get update
ADD clean_fasta_wrapper.sh /usr/bin
RUN chmod +x /usr/bin/clean_fasta_wrapper.sh

ENTRYPOINT ["/usr/bin/clean_fasta_wrapper.sh"]
