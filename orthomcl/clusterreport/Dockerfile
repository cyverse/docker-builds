## Dockerfile
FROM ubuntu:16.04
MAINTAINER Amanda Cooksey	
LABEL Description="A utility to generate species-centric files and a report on the number clusters produced by OrthoMCL"

# Install all the updates and download dependencies
RUN apt-get update && apt-get install -y git

ADD clusterReport.pl /usr/bin

# Change the permissions and the path for the script
RUN chmod +x /usr/bin/clusterReport.pl

# Entrypoint
ENTRYPOINT ["/usr/bin/clusterReport.pl"]

RUN mkdir /work-dir

WORKDIR /work-dir

