## Dockerfile
FROM ubuntu:18.04
MAINTAINER Amanda Cooksey	
LABEL Description="rRNAFilter"

# Install all the updates and download dependencies
RUN apt-get update && \
    apt-get install -y \
    openjdk-11-jre \ 
    git \
    wget \
    bzip2 \
    zip \
    nano

WORKDIR /usr/local/bin


ADD rRNAFilter/count \
    rRNAFilter/java-getopt-1.0.14.jar \
    rRNAFilter/rRNAFilter_commandline.jar\
    rRNAFilter/kanalyze.jar \
    wrapper.sh \
    /usr/local/bin/

RUN chmod a+x /usr/local/bin/wrapper.sh

#RUN mkdir /workdir

#WORKDIR /workdir

#ENV PATH=$PATH:/workdir

ENTRYPOINT ["bash", "/usr/local/bin/wrapper.sh"]
