FROM ubuntu:14.04.3
MAINTAINER Upendra Devisetty <upendra@cyverse.org>
LABEL Description "This Dockerfile is for modifying the coordinates of gff file"

RUN apt-get update
ADD Modify_gff_coordinates.sh /usr/bin
RUN chmod +x /usr/bin/Modify_gff_coordinates.sh

ENTRYPOINT ["/usr/bin/Modify_gff_coordinates.sh"]
