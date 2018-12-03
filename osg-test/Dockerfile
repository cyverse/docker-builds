FROM ubuntu:xenial

RUN mkdir /cvmfs /work
WORKDIR /work

# Install some prerequisites.
RUN apt-get update \
    && apt-get install -y jq lsb wget curl apt-transport-https python2.7 python-requests libfuse2

# Define the iRODS package.
ENV ICMD_BASE="https://files.renci.org/pub/irods/releases/4.1.10/ubuntu14"
ENV ICMD_PKG="irods-icommands-4.1.10-ubuntu14-x86_64.deb"

# Install icommands.
RUN curl -o "$ICMD_PKG" "$ICMD_BASE/$ICMD_PKG" \
    && dpkg -i "$ICMD_PKG" \
    && rm -f "$ICMD_PKG"

# Install the wrapper script and the script to upload the output files.
COPY wrapper /usr/bin/wrapper
COPY send-job-status-update /usr/bin/send-job-status-update
COPY ticketParser.sh /usr/bin/ticketParser.sh
COPY evalTicket.sh /usr/bin/evalTicket.sh

# Make the wrapper script the default command.
CMD ["wrapper"]
