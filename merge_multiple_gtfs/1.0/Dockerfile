FROM ubuntu:16.04
MAINTAINER Upendra Devsietty <upendra@cyverse.org>
LABEL Description "This Dockerfile is for merge_multiple_gtfs-1.0"

RUN apt-get update && apt-get install -y g++ \
					 wget \
					 make \
					 zlib1g-dev \
					 python
# Cufflinks
RUN wget -O- http://cole-trapnell-lab.github.io/cufflinks/assets/downloads/cufflinks-2.2.1.Linux_x86_64.tar.gz | tar xzvf -

# bedtools
RUN wget -O- https://github.com/arq5x/bedtools2/releases/download/v2.25.0/bedtools-2.25.0.tar.gz | tar xzvf -
WORKDIR bedtools2/
RUN make

# Wrapper script
ADD merge_multiple_gtfs.sh /usr/bin
RUN chmod +x /usr/bin/merge_multiple_gtfs.sh

ENV PATH /bedtools2/bin/:$PATH
ENV PATH /cufflinks-2.2.1.Linux_x86_64/:$PATH

ENTRYPOINT ["merge_multiple_gtfs.sh"]
