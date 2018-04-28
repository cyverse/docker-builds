FROM r-base:3.4.4
MAINTAINER Upendra Devisetty <upendra@cyverse.org>
LABEL Description "This Dockerfile is for merge_htseq:1.0"

RUN apt-get update

RUN Rscript -e 'install.packages("getopt", dependencies = TRUE);'

Add merge_htseq.R /usr/bin

RUN chmod +x /usr/bin/merge_htseq.R

ENTRYPOINT ["merge_htseq.R"]

CMD ["merge_htseq.R", "-h"]