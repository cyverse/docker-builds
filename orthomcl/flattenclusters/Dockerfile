## Dockerfile
FROM ubuntu:16.04
MAINTAINER Amanda Cooksey	
LABEL Description="A utility to map original Fasta Headers to clusterReport and/or queryOrthoMCL output."

# Install all the updates and download dependencies
RUN apt-get update && apt-get install -y git

ADD flattenClusters.pl /usr/bin

# Change the permissions and the path for the script
RUN chmod +x /usr/bin/flattenClusters.pl

# Entrypoint
ENTRYPOINT ["/usr/bin/flattenClusters.pl"]

RUN mkdir /work-dir

WORKDIR /work-dir

