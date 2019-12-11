## Dockerfile
FROM ubuntu:16.04
MAINTAINER Amanda Cooksey	
LABEL Description="A utility to add unclustered sequences to OrthoMCL output."

# Install all the updates and download dependencies
RUN apt-get update && apt-get install -y git

ADD appendUnclustered.pl /usr/bin

# Change the permissions and the path for the script
RUN chmod +x /usr/bin/appendUnclustered.pl

# Entrypoint
ENTRYPOINT ["/usr/bin/appendUnclustered.pl"]

RUN mkdir /work-dir

WORKDIR /work-dir

