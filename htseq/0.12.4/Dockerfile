# Set the base image to Ubuntu
FROM ubuntu:18.04

# File Author / Maintainer
MAINTAINER Amanda Cooksey <amcooksey@cyverse.org>

# Update the repository sources list
RUN apt-get update

RUN DEBIAN_FRONTEND="noninteractive" apt-get -y install tzdata

# Install Latex
RUN apt-get install --yes build-essential python-dev python-numpy python-matplotlib python-pip zlib1g-dev wget

# Download and uncompress HTSeq archive
WORKDIR /tmp
RUN wget https://files.pythonhosted.org/packages/c4/04/b9b0c5514dcd09e64481e8ebc242aef162646b6de956ffb44595d1de0f69/HTSeq-0.12.4.tar.gz
RUN tar xvf HTSeq-0.12.4.tar.gz
WORKDIR HTSeq-0.12.4

# Install HTSeq
RUN python setup.py build && python setup.py install

# Install pysam
#RUN pip install pysam==0.9.0

# Cleanup
RUN apt-get clean
#RUN rm -rf /tmp/*

# Wrapper script
ADD htseqcount_wrapper.sh /usr/bin
RUN chmod +x /usr/bin/htseqcount_wrapper.sh

# ENTRYPOINT
ENTRYPOINT ["htseqcount_wrapper.sh"]
CMD ["-h"]
