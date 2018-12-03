FROM ubuntu:xenial

RUN mkdir /cvmfs /work
WORKDIR /work

# Install some prerequisites.
RUN apt-get update \
    && apt-get install -y lsb curl wget apt-transport-https python2.7 python-requests libfuse2

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

# Make the wrapper script the default command.
CMD ["wrapper"]
