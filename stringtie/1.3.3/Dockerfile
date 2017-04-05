# StringTie 
# VERSION             1.3.3
# 

FROM      ubuntu:14.04.3 
MAINTAINER Kapeel Chougule <kchougul@cshl.edu>
LABEL Description="This image is used for running StringTie-1.3.3 assembly and quantification tool for RNA seq data"

RUN apt-get update && apt-get install -y build-essential wget python rsync
ADD StringTie_pl.pl /usr/bin/ 
ADD prepDE.py /usr/bin/prepDE.py 
RUN [ "chmod", "+x", "/usr/bin/StringTie_pl.pl" ] 
RUN [ "chmod", "+x", "/usr/bin/prepDE.py" ] 
RUN wget http://ccb.jhu.edu/software/stringtie/dl/stringtie-1.3.3.Linux_x86_64.tar.gz \
	&& tar -xvf stringtie-1.3.3.Linux_x86_64.tar.gz \
	&& mv stringtie-1.3.3.Linux_x86_64 stringtie-1.3.3 \
	&& cp stringtie-1.3.3/stringtie /usr/bin/

ENTRYPOINT ["StringTie_pl.pl"]


