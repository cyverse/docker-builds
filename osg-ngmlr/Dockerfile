FROM ubuntu:xenial
MAINTAINER Amanda Cooksey <amcooksey@email.arizona.edu>

RUN mkdir /cvmfs /work

RUN apt-get update \
    && apt-get install -y wget gcc cmake libz-dev curl lsb apt-transport-https python3 python-requests libfuse2 make libpcre3-dev


# Install ngmlr
RUN wget https://github.com/philres/ngmlr/releases/download/v0.2.7/ngmlr-0.2.7-linux-x86_64.tar.gz
RUN tar xvzf ngmlr-0.2.7-linux-x86_64.tar.gz
RUN mv ngmlr-0.2.7/ngmlr /usr/bin/

WORKDIR /work

# Define the iRODS package.
ENV ICMD_BASE="https://files.renci.org/pub/irods/releases/4.1.10/ubuntu14"
ENV ICMD_PKG="irods-icommands-4.1.10-ubuntu14-x86_64.deb"

# Install icommands.
RUN curl -o "$ICMD_PKG" "$ICMD_BASE/$ICMD_PKG" \
        && dpkg -i "$ICMD_PKG" \
        && rm -f "$ICMD_PKG"

# Install the wrapper script and the script to upload the output files.
ADD wrapper /usr/bin/wrapper
ADD upload-files /usr/bin/upload-files

RUN chmod +x /usr/bin/upload-files /usr/bin/wrapper

# Make the wrapper script the default command.
CMD ["wrapper"]
