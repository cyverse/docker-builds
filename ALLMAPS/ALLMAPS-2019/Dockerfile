## Dockerfile
FROM ubuntu:18.04
MAINTAINER Amanda Cooksey	
LABEL Description="ALLMAPS--uses UCSC LiftOver--not permitted for commericial use."

ENV DEBIAN_FRONTEND noninteractive

# Install all the updates and download dependencies
RUN apt-get update
RUN apt-get install --yes  \
	apt-utils \
	python \
	python-pip \
	curl \
	wget \
	nano \
	libfontconfig1 \
	font-manager \
	git \
	texlive \
	texlive-latex-extra \
	texlive-latex-recommended

RUN pip install --user biopython==1.70 numpy deap networkx matplotlib jcvi 

RUN wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/liftOver

RUN wget --no-check-certificate https://www.dropbox.com/s/onsjieazu0uytgk/ALLMAPS-install.sh && \
	sh ALLMAPS-install.sh && \
	cp concorde faSize liftOver /usr/bin/

RUN mkdir /usr/bin/allmaps

WORKDIR /usr/bin/allmaps


##THERE IS A BUG THAT WAS FIXED? IN SEPT--BE SURE TO GET LATEST VERSION--THAT IS MY CURRENT ISSUE
RUN git clone git://github.com/tanghaibao/jcvi.git

ENV PYTHONPATH=/usr/bin/allmaps:$PYTHONPATH

ENV PATH=/root/.local/lib/python2.7/site-packages:$PATH



# Entrypoint
ENTRYPOINT ["python","-m","jcvi.assembly.allmaps"]

RUN mkdir /work-dir

WORKDIR /work-dir
