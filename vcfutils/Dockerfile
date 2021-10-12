## Dockerfile
FROM ubuntu:20.04
MAINTAINER Amanda Cooksey	
LABEL Description="varfilter.py"

# Install all the updates and download dependencies
RUN apt-get update && \
    apt-get install -y \
    git \
    wget \
    bzip2 

RUN echo 'export PATH=/opt/conda/bin:$PATH' > /etc/profile.d/conda.sh && \
    wget --quiet https://repo.continuum.io/miniconda/Miniconda2-4.0.5-Linux-x86_64.sh -O ~/miniconda.sh && \
    /bin/bash ~/miniconda.sh -b -p /opt/conda && \
    rm ~/miniconda.sh


# give write permissions to conda folder
RUN chmod 777 -R /opt/conda/

ENV PATH=$PATH:/opt/conda/bin

RUN conda config --add channels bioconda

RUN conda upgrade conda

RUN conda install -c conda-forge -c bioconda bcftools==1.13


# add blastp_uniprot--i added this
ADD vcfutils.pl /usr/bin

# Change the permissions and the path for the wrapper script
RUN chmod +x /usr/bin/vcfutils.pl

# Entrypoint
ENTRYPOINT ["/usr/bin/vcfutils.pl"]

