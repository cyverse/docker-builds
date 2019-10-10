# STAR 
# VERSION              2.5.3.a
#

FROM      ubuntu:14.04.3 
MAINTAINER Kapeel Chougule <kchougul@cshl.edu>
LABEL Description="This image is used for running STAR RNA seq alignemnt tool version 2.5.3.a"

RUN apt-get update && apt-get install -y build-essential wget git curl unzip zlib1g-dev
ADD STAR_pl.pl /usr/bin/ 
RUN [ "chmod", "+x", "/usr/bin/STAR_pl.pl" ]
RUN wget https://github.com/alexdobin/STAR/archive/2.5.3a.tar.gz
RUN tar zxvf 2.5.3a.tar.gz
RUN cp /STAR-2.5.3a/bin/Linux_x86_64/STAR /usr/bin/

ENTRYPOINT ["STAR_pl.pl"]
