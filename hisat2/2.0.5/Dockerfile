FROM ubuntu:14.04.3
MAINTAINER Kapeel Chougule <kchougul@cshl.edu>
LABEL Description "This Dockerfile is used for hisat2 tool"
RUN apt-get update && apt-get install -y \
   build-essential \
   git \
   python \
   samtools
ENV BINPATH /usr/bin
ENV SRCPATH /usr/src
ENV HISAT2GIT https://github.com/infphilo/hisat2.git
ENV HISAT2PATH $SRCPATH/hisat2
RUN mkdir -p $SRCPATH
WORKDIR $SRCPATH

RUN git clone "$HISAT2GIT" \
   && cd $HISAT2PATH \
RUN  make -C $HISAT2PATH \
   && cp $HISAT2PATH/hisat2 $BINPATH \
   && cp $HISAT2PATH/hisat2-* $BINPATH

ADD Hisat2_align.pl $BINPATH
RUN chmod +x $BINPATH/Hisat2_align.pl

ENTRYPOINT ["/usr/bin/Hisat2_align.pl"]
