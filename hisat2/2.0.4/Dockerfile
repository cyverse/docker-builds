FROM ubuntu:14.04.3
MAINTAINER Kapeel Chougule
RUN apt-get update && apt-get install -y \
   build-essential \
   git \
   python
ENV BINPATH /usr/bin
ENV SRCPATH /usr/src
ENV HISAT2GIT https://github.com/infphilo/hisat2.git
ENV HISAT2PATH $SRCPATH/hisat2
RUN mkdir -p $SRCPATH
ADD HISAT2_align.pl $BINPATH
WORKDIR $SRCPATH
RUN git clone "$HISAT2GIT" \
   && cd $HISAT2PATH \
RUN  make -C $HISAT2PATH \
   && cp $HISAT2PATH/hisat2 $BINPATH \
   && cp $HISAT2PATH/hisat2-* $BINPATH
ENTRYPOINT ["HISAT2_align.pl"]
