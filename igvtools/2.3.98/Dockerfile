FROM java:7
MAINTAINER Upendra Devisetty <upendra@cyverse.org>
LABEL Description "This Dockerfile is used for igvtools 2.3.98"

RUN apt-get update && apt-get install -y wget

WORKDIR /igv

RUN wget http://data.broadinstitute.org/igv/projects/downloads/2.3/igvtools_2.3.98.zip
RUN unzip igvtools_2.3.98.zip 

WORKDIR /

ENTRYPOINT ["java", "-jar", "/igv/IGVTools/igvtools.jar"]
