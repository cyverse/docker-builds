FROM ubuntu:16.04
MAINTAINER Upendra Devisetty <upedra@cyverse.org>
LABEL Description "This Dockerfile is for Synapse Client 1.6.1"
ENV PACKAGES python-dev git python-setuptools python-pip

ENV BRANCH=develop
ENV VERSION=6ba6a3ebde81fe8ed4d0c231ab42c613aa03334f

RUN apt-get update && \
    apt-get install -y --no-install-recommends ${PACKAGES}

RUN git clone -b ${BRANCH} git://github.com/Sage-Bionetworks/synapsePythonClient.git && \
    cd synapsePythonClient && \
    git checkout ${VERSION} && \
    python setup.py develop

ADD synapse_wrapper.sh /usr/bin
RUN chmod +x /usr/bin/synapse_wrapper.sh

ENTRYPOINT ["synapse_wrapper.sh"]
CMD ["-h"]
