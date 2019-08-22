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
        && ./configure LIBnn=lib --with-x --with-readline --prefix=/opt/R3.6 --enable-R-shlib \
            --with-cairo --with-libpng --with-blas --with-lapack \
        && make \
        && make install \
        && ln -s /opt/R3.6/bin/r /opt/conda/bin/R \
        && ln -s /opt/R3.6/bin/Rscript /opt/conda/bin/Rscript

# Place the new Rprofile.site file.
COPY --chown=jovyan:users Rprofile.site /opt/R3.6/lib/R/etc/Rprofile.site

# Install common R libraries.
RUN /opt/R3.6/bin/Rscript -e "install.packages('openssl')" \
        && /opt/R3.6/bin/Rscript -e "install.packages('curl')" \
        && /opt/R3.6/bin/Rscript -e "install.packages('httr')" \
        && /opt/R3.6/bin/Rscript -e "install.packages('ggplot2')" \
        && /opt/R3.6/bin/Rscript -e "install.packages('BiocManager')" \
        && /opt/R3.6/bin/Rscript -e "install.packages('IRkernel')"

# Restore the environment after R compilation
USER root
ENV PATH=$OLD_PATH
RUN chmod 755 /opt && \
rm -r /opt/r-tmp

#############################################
##  $NB_USER                               ##
##      Update the theme                   ##
#############################################

USER $NB_USER

# Copy the logo and favicon
COPY logo.png favicon.ico /opt/conda/lib/python3.6/site-packages/notebook/static/base/images/

# Install the custom theme
COPY custom.css custom.js theme.css background.jpg /opt/conda/lib/python3.6/site-packages/notebook/static/custom/

#############################################
##  ROOT                                   ##
##      Configure nbextensions             ##
#############################################

USER root

# Install GenePattern nbextensions in default Python environment
RUN pip install -v nbtools genepattern-notebook jupyter_wysiwyg

# Enable the extensions
RUN jupyter nbextension enable --sys-prefix --py nbtools \
        && jupyter nbextension enable --sys-prefix --py genepattern \
        && jupyter nbextension enable --sys-prefix --py jupyter_wysiwyg

# Enable Notebook Repository frontend
RUN git clone https://github.com/genepattern/notebook-repository-frontend.git \
        && cd notebook-repository-frontend \
        && jupyter nbextension install repo \
        && jupyter nbextension enable --sys-prefix repo/js/main \
        && jupyter nbextension enable --sys-prefix --section=tree repo/js/main \
        && rm -r ../notebook-repository-frontend

# Enable Notebook Repository hints
RUN git clone https://github.com/genepattern/notebook-tour.git \
        && cd notebook-tour \
        && jupyter nbextension install hints \
        && jupyter nbextension enable --sys-prefix --section=tree hints/js/main \
        && rm -r ../notebook-tour

# Install the collapsible_headings extension
ENV JC_BASE=https://github.com/ipython-contrib/jupyter_contrib_nbextensions/raw/master/src
ENV CH_BASE=$JC_BASE/jupyter_contrib_nbextensions/nbextensions/collapsible_headings
ENV CH_DEST=/opt/conda/share/jupyter/nbextensions/collapsible_headings
ADD $CH_BASE/main.js $CH_DEST
ADD $CH_BASE/main.css $CH_DEST
COPY collapsible_headings.json $CH_DEST
