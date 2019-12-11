## Dockerfile
FROM ubuntu:16.04
MAINTAINER Amanda Cooksey	
LABEL Description="A utility to generate custom homolog sets based on the copy number of genes for each species."

# Install all the updates and download dependencies
RUN apt-get update && apt-get install -y git

ADD queryOrthoMCL.pl /usr/bin

# Change the permissions and the path for the script
RUN chmod +x /usr/bin/queryOrthoMCL.pl

# Entrypoint
ENTRYPOINT ["/usr/bin/queryOrthoMCL.pl"]

RUN mkdir /work-dir

WORKDIR /work-dir

