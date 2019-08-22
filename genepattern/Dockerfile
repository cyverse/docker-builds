# This Dockerfile was copied and modified from the original GenePattern Notebook Dockerfile:
# https://raw.githubusercontent.com/genepattern/notebook-docker/master/genepattern-notebook/Dockerfile

FROM cyversevice/jupyterlab-scipy:latest

USER root

# Install prerequisites for compiling R
RUN echo "deb https://cloud.r-project.org/bin/linux/ubuntu bionic-cran35/" >> /etc/apt/sources.list \
        && apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9 \
        && apt-get update \
        && apt-get install -y r-base-dev libcairo2-dev libxt-dev xvfb libcurl4-openssl-dev

# Set up the environment for R compilation (take conda off the PATH, it screws up compiling R)
ENV OLD_PATH=$PATH
ENV PATH=/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin
RUN chmod 777 /opt
USER $NB_USER

# Download and compile R.
ADD --chown=jovyan:users http://cran.r-project.org/src/base/R-3/R-3.6.1.tar.gz /opt/r-tmp/
RUN cd /opt/r-tmp \
        && tar xvf R-3.6.1.tar.gz \
        && cd R-3.6.1 \
        && ./configure LIBnn=lib --with-x --with-readline --prefix=/opt/r3.6 --enable-R-shlib \
            --with-cairo --with-libpng --with-blas --with-lapack \
        && make \
        && make install \
        && ln -s /opt/R3.6/bin/r /opt/conda/bin/R \
        && ln -s /opt/R3.6/bin/Rscript /opt/conda/bin/Rscript
