## Dockerfile
FROM ubuntu:16.04
MAINTAINER Amanda Cooksey	
LABEL Description="Dockerized version of fastaRename for the OrthoMCL pipeline"

# Install all the updates and download dependencies
RUN apt-get update && apt-get install -y git

ADD fastaRename.pl /usr/bin

# Change the permissions and the path for the script
RUN chmod +x /usr/bin/fastaRename.pl

# Entrypoint
ENTRYPOINT ["/usr/bin/fastaRename.pl"]

RUN mkdir /work-dir

WORKDIR /work-dir

