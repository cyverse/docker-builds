FROM ubuntu:16.04
MAINTAINER Upendra Devisetty <upendra@cyverse.org>
LABEL "version":1.0 "Description":"This Dockerfile is for SNP Discovery"

RUN apt-get update && apt-get install -y cpanminus

COPY SNP_Discovery.pl /
RUN chmod +x /SNP_Discovery.pl && cp /SNP_Discovery.pl /usr/bin

ENTRYPOINT ["SNP_Discovery.pl"]
CMD ["-h"] 
