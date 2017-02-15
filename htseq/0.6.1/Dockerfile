# Set the base image to Ubuntu
FROM ubuntu:14.04.3

# File Author / Maintainer
MAINTAINER Upendra Devisetty <upendra@cyverse.org>

# Update the repository sources list
RUN apt-get update

# Install Latex
RUN apt-get install --yes build-essential python2.7-dev python-numpy python-matplotlib python-pip zlib1g-dev wget

# Download and uncompress HTSeq archive
WORKDIR /tmp
RUN wget https://pypi.python.org/packages/72/0f/566afae6c149762af301a19686cd5fd1876deb2b48d09546dbd5caebbb78/HTSeq-0.6.1.tar.gz
RUN tar xvf HTSeq-0.6.1.tar.gz
WORKDIR HTSeq-0.6.1

# Install HTSeq
RUN python setup.py build && python setup.py install

# Install pysam
RUN pip install 'pysam'

# Cleanup
RUN apt-get clean
RUN rm -rf /tmp/*

# Wrapper script
ADD htseqcount_wrapper.sh /usr/bin
RUN chmod +x /usr/bin/htseqcount_wrapper.sh

# ENTRYPOINT
ENTRYPOINT ["htseqcount_wrapper.sh"]
CMD ["-h"]
