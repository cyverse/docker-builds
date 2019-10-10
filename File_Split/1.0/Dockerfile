FROM ubuntu:16.04
MAINTAINER Upendra Devisetty <upendra@cyverse.org>
LABEL "version":"1.0" "Description":"This Dockerfile is for file split"

RUN apt-get update
COPY file_split.sh /usr/bin
RUN chmod +x /usr/bin/file_split.sh

ENTRYPOINT ["file_split.sh"]
CMD ["-h"]
