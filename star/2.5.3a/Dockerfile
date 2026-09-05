FROM ubuntu:22.04

LABEL maintainer="Kapeel Chougule <kchougul@cshl.edu>"
LABEL description="This image is used for running the STAR RNA-seq alignment tool version 2.5.3a"

ARG STAR_VERSION=2.5.3a

# perl runs the wrapper script; gzip feeds STAR's --readFilesCommand;
# libgomp1 and zlib1g are STAR's runtime dependencies.
RUN apt-get update \
    && apt-get install -y --no-install-recommends \
        ca-certificates \
        gzip \
        libgomp1 \
        perl \
        wget \
        zlib1g \
    && rm -rf /var/lib/apt/lists/*

# Extract only the precompiled binary, so the source tree does not
# end up as a layer in the image.
RUN wget -q -O /tmp/STAR.tar.gz \
        "https://github.com/alexdobin/STAR/archive/refs/tags/${STAR_VERSION}.tar.gz" \
    && tar -xzf /tmp/STAR.tar.gz -C /tmp \
        "STAR-${STAR_VERSION}/bin/Linux_x86_64/STAR" \
    && install -m 0755 "/tmp/STAR-${STAR_VERSION}/bin/Linux_x86_64/STAR" /usr/local/bin/STAR \
    && rm -rf /tmp/STAR.tar.gz "/tmp/STAR-${STAR_VERSION}"

COPY STAR_pl.pl /usr/local/bin/STAR_pl.pl
RUN chmod +x /usr/local/bin/STAR_pl.pl

WORKDIR /work-dir

ENTRYPOINT ["STAR_pl.pl"]

# Build
# docker build -t cyverse/star:2.5.3a .
# Test run
# docker run --rm -v $(pwd):/work-dir -w /work-dir cyverse/star:2.5.3a \
#     --file_query reads_1.fastq.gz --user_database genome.fasta \
#     --user_annotation genes.gtf --file_type SE
