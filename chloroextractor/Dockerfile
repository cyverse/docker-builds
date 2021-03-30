ARG osversion=v1.0.0
FROM chloroextractorteam/chloroextractor-dockerbase:${osversion}

ARG VERSION=v1.0.9
ARG VCS_REF
ARG BUILD_DATE
ARG osversion

RUN echo "osversion: "${osversion}", VCS_REF: "${VCS_REF}", BUILD_DATE: "${BUILD_DATE}", VERSION: "${VERSION}

LABEL maintainer="frank.foerster@ime.fraunhofer.de" \
      description="Dockerfile providing the image for our chloroExtractor container" \
      version=${VERSION} \
      org.label-schema.vcs-ref=${VCS_REF} \
      org.label-schema.build-date=${BUILD_DATE} \
      org.label-schema.vcs-url="https://github.com/chloroExtractorTeam/chloroExtractor.git"

RUN git clone --recursive \
              --branch ${VERSION} \
        https://github.com/chloroExtractorTeam/chloroExtractor.git \
        /opt/chloroExtractor/ && \
    rm -rf /opt/chloroExtractor/.git
ENV PATH "/opt/chloroExtractor/bin/:$PATH"

#VOLUME /data

#WORKDIR /data

ENTRYPOINT ["/opt/chloroExtractor/bin/ptx"]

CMD ["--help"]
