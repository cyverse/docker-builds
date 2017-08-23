FROM debian:wheezy
MAINTAINER Upendra Devisetty <upendra@cyverse.org>
LABEL Description "This Dockerfile is for blastp-2.6.0"

RUN apt-get update && apt-get install -y --no-install-recommends \
		curl \
		python-minimal=2.7.3-4+deb7u1 \
	&& rm -rf /var/lib/apt/lists/*

# Download & install BLAST
RUN mkdir /opt/blast \
      && curl ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ncbi-blast-2.6.0+-x64-linux.tar.gz \
      | tar -zxC /opt/blast --strip-components=1
ENV PATH /opt/blast/bin:$PATH

# Add wrapper script
ADD blastp_wrapper.sh /usr/bin
RUN chmod +x /usr/bin/blastp_wrapper.sh

# Entrypoint
ENTRYPOINT ["blastp_wrapper.sh"]
