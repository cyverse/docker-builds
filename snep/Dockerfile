FROM ubuntu:14.04

MAINTAINER Andre Mercer <amercer@iplantcollaborative.org>

# Run updates
RUN apt-get update
RUN apt-get upgrade -y

COPY SNeP1.1 /usr/local/bin/

ENTRYPOINT ["SNeP1.1"]
