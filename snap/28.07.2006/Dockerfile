FROM ubuntu:14.04.3
MAINTAINER Upendra Devisetty <upendra@cyverse.org>
LABEL Description "This Dockerfile is used for SNAP gene prediction tool version version 2006-07-28"

# Run updates
RUN apt-get update
RUN apt-get -y install wget make gcc

# SNAP
RUN wget http://korflab.ucdavis.edu/Software/snap-2013-11-29.tar.gz
RUN tar xvf snap-2013-11-29.tar.gz
WORKDIR /snap
RUN make
RUN chmod +x snap && cp snap /usr/bin

# Add wrapper script
ADD snap_28_07_gene_prediction.sh /usr/bin
RUN chmod +x /usr/bin/snap_28_07_gene_prediction.sh

# Entrypoint
ENTRYPOINT ["/usr/bin/snap_28_07_gene_prediction.sh"]