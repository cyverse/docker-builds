## Dockerfile
FROM ubuntu:16.04
MAINTAINER Amanda Cooksey	
LABEL Description="A Tool for Statistical Analysis and Visualization of High-Throughput Biological Data Using Gene Ontology. Additional scripts are used to ensure Ontologizer, the GO terminology and structure file (.obo) and the GAF file are up to date."

# Install all the updates and download dependencies
RUN apt-get update && apt-get install -y \
	wget \
	tar \
	nano \
	default-jre-headless

RUN mkdir /work-dir /ontologizer

ADD 2020_Ontologizer.sh Ontologizer.jar /ontologizer/

ENV PATH=$PATH:/ontologizer

# Change the permissions and the path for the script
RUN chmod +x /ontologizer/2020_Ontologizer.sh /ontologizer/Ontologizer.jar

# Entrypoint
ENTRYPOINT ["/ontologizer/2020_Ontologizer.sh"]

WORKDIR /work-dir

