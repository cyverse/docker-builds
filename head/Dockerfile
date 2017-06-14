FROM alpine
MAINTAINER Dennis Roberts <dennis@cyverse.org>

RUN apk add --update "perl=5.24.0-r0" \
    && rm -rf /var/cache/apk/*

COPY headWrapper.pl /usr/local/bin

ENTRYPOINT ["headWrapper.pl"]
