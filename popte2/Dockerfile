From ubuntu:18.04
MAINTAINER Upendra Devisetty <upendra@cyverse.org>
LABEL version="v1.10.04" description="This tool is used for PoPoolation TE2"

RUN apt-get update && apt-get install -y \
			wget \
			openjdk-8-jdk \

RUN echo 'export PATH=/opt/conda/bin:$PATH' > /etc/profile.d/conda.sh && \
    wget --quiet https://repo.continuum.io/miniconda/Miniconda2-4.3.14-Linux-x86_64.sh -O ~/miniconda.sh && \
    /bin/bash ~/miniconda.sh -b -p /opt/conda && \
    rm ~/miniconda.sh

RUN conda install -y -c bioconda bwa=0.7.17
RUN conda install -y -c bioconda samtools=1.9
ENV PATH /opt/conda/bin:$PATH

COPY run-map-subs-popte-wrapper.sh popte2-v1.10.04.jar /

RUN chmod +x /run-map-subs-popte-wrapper.sh 

ENTRYPOINT ["run-map-subs-popte-wrapper.sh"]
CMD ["-h"]

