# StringTie 
# VERSION             1.2.4
# 

FROM      ubuntu:14.04.3 
MAINTAINER Kapeel Chougule <kchougul@cshl.edu>
LABEL Description="This image is used for running StringTie-1.2.4 assembly and quantification tool for RNA seq data"

RUN apt-get update && apt-get install -y build-essential wget python 
ADD StringTie.pl /usr/bin/ 
ADD prepDE.py /usr/bin/prepDE.py 
RUN [ "chmod", "+x", "/usr/bin/StringTie.pl" ] 
RUN [ "chmod", "+x", "/usr/bin/prepDE.py" ] 
RUN wget http://ccb.jhu.edu/software/stringtie/dl/stringtie-1.2.4.Linux_x86_64.tar.gz \
	&& tar -xvf stringtie-1.2.4.Linux_x86_64.tar.gz \
	&& mv stringtie-1.2.4.Linux_x86_64 stringtie-1.2.4 \
	&& cp stringtie-1.2.4/stringtie /usr/bin/

ENTRYPOINT ["StringTie.pl"]


