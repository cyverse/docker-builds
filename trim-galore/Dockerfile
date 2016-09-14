FROM ubuntu:14.04

MAINTAINER Andre Mercer <amercer@iplantcollaborative.org>

# Run updates
RUN apt-get clean
RUN apt-get update
RUN apt-get upgrade -y
RUN apt-get install wget -y
RUN apt-get install unzip -y

RUN wget http://www.bioinformatics.babraham.ac.uk/projects/trim_galore/trim_galore_v0.4.1.zip
RUN unzip trim_galore_v0.4.1.zip

ENV PATH /trim_galore_zip:$PATH

ENTRYPOINT ["trim_galore"]
