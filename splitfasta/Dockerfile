## Dockerfile
FROM ubuntu:16.04
MAINTAINER Amanda Cooksey	
LABEL Description="splitfasta.pl"

# Install all the updates and download dependencies
RUN apt-get update && apt-get install -y git

ADD splitfasta.pl /usr/bin

# Change the permissions and the path for the script
RUN chmod +x /usr/bin/splitfasta.pl

# Entrypoint
ENTRYPOINT ["/usr/bin/splitfasta.pl"]

RUN mkdir /work-dir

WORKDIR  /work-dir

