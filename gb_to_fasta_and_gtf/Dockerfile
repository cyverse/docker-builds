FROM ubuntu:latest

MAINTAINER Michael Thon <mike.thon@gmail.com>

RUN apt-get update
RUN apt-get install -y emboss python
ADD gb_to_fasta_and_gtf.py .
ENTRYPOINT gb_to_fasta_and_gtf.py
