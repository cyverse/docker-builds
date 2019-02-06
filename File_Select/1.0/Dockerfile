FROM ubuntu:16.04
MAINTAINER Upendra Devisetty <upendra@cyverse.org>
LABEL "version":"1.0" "Description":"This Dockerfile is for file select"

RUN apt-get update
COPY file_select.sh /usr/bin
RUN chmod +x /usr/bin/file_select.sh

ENTRYPOINT ["file_select.sh"]
CMD ["-h"]
