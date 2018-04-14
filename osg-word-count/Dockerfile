FROM ubuntu:xenial

RUN mkdir /cvmfs /work
WORKDIR /work

# Install some prerequisites.
RUN apt-get update \
    && apt-get install -y lsb wget apt-transport-https python2.7

# Install icommands.
RUN wget -qO - https://packages.irods.org/irods-signing-key.asc | apt-key add - \
    && echo "deb [arch=amd64] https://packages.irods.org/apt/ xenial main" > /etc/apt/sources.list.d/renci-irods.list \
    && apt-get update \
    && apt-get install -y irods-icommands
