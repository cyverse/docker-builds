FROM ubuntu:14.04.3
MAINTAINER Upendra Devisetty <upendra@cyverse>
LABEL Description "This Dockerfile is for SPADES-3.9.0 version"

RUN apt-get update && apt-get install -y wget python

ADD http://spades.bioinf.spbau.ru/release3.9.0/SPAdes-3.9.0-Linux.tar.gz /tmp/
RUN tar zxf /tmp/SPAdes-3.9.0-Linux.tar.gz

WORKDIR /SPAdes-3.9.0-Linux
RUN mv bin/* /usr/local/bin/
RUN mv share/* /usr/local/share
WORKDIR /

ENTRYPOINT ["/usr/local/bin/spades.py"]
CMD ["/usr/local/bin/spades.py -h"]

# Build
# sudo docker build -t"=ubuntu/spades:3.9.0" .
# Run without arguments
# sudo docker run ubuntu/spades:3.9.0 -h
