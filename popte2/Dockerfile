From ubuntu:18.04
MAINTAINER Upendra Devisetty <upendra@cyverse.org>
LABEL version="v1.10.04" description="This tool is used for PoPoolation TE2"

RUN apt-get update && apt-get install -y \
                        wget \
                        openjdk-8-jdk

RUN mkdir -p /tmp/conda-build && wget -nv https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh && bash Miniconda3-latest-Linux-x86_64.sh -b -p /miniconda3
RUN cp /miniconda3/bin/conda /miniconda3/bin/conda3
ENV PATH=/miniconda3/bin:$PATH
RUN conda3 update conda
RUN conda3 config --add channels bioconda
RUN conda3 install samtools -c bioconda -y
RUN conda3 install bwa -c bioconda -y

COPY run-map-subs-popte-wrapper.sh popte2-v1.10.04.jar /
RUN chmod +x /run-map-subs-popte-wrapper.sh

ENTRYPOINT ["/run-map-subs-popte-wrapper.sh"]
CMD ["-h"]
