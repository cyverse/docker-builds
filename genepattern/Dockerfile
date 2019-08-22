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
ADD $CH_BASE/main.js $CH_DEST/
ADD $CH_BASE/main.css $CH_DEST/
COPY collapsible_headings.json $CH_DEST/

# Install the toc2 extension
ENV TOC2_BASE=$JC_BASE/jupyter_contrib_nbextensions/nbextensions/toc2
ENV TOC2_DEST=/opt/conda/share/jupyter/nbextensions/toc2
ADD $TOC2_BASE/main.js $TOC2_DEST/
ADD $TOC2_BASE/main.css $TOC2_DEST/
ADD $TOC2_BASE/toc2.js $TOC2_DEST/
COPY toc2.json $TOC2_DEST/
RUN sed -i "s/number_sections: true/number_sections: false/" $TOC2_DEST/toc2.js \
        && sed -i "s/moveMenuLeft: true/moveMenuLeft: false/" $TOC2_DEST/toc2.js

# Install the igv.js extension
RUN pip install -v igv-jupyter==0.9.7 \
        && jupyter serverextension enable --sys-prefix --py igv \
        && jupyter nbextension install --sys-prefix --py igv \
        && jupyter nbextension enable --sys-prefix --py igv

# Install the cytoscape extension
RUN pip install -v cyjupyter==0.2.0 \
        && jupyter nbextension enable --sys-prefix --py cyjupyter

##########################################
##  $NB_USER                            ##
##      Create the Python environments  ##
##########################################

USER $NB_USER

# Create the Python 3.6 environment
RUN conda create -y --name python3.6 python=3.6 anaconda ipykernel==5.1.0

# Install the GenePattern Notebook Environment packages. Running these commands in a new shell isolates
# the conda environment.
RUN /bin/bash -c "source activate python3.6 && \
    conda install -y -c conda-forge 'biopython==1.70' 'igraph==0.7.1' 'babel==2.4.0' 'louvain==0.6.1'"
RUN /bin/bash -c "source activate python3.6 && \
    pip install nbtools genepattern-python genepattern-notebook jupyter_wysiwyg 'igv-jupyter==0.9.7' \
    'cyjupyter==0.2.0' 'ccalnoir==2.6' 'cuzcatlan==0.9.0' 'ndex2==1.2.0.*' 'anndata==0.5.10' 'enum34==1.1.6' \
    'ijson==2.3' 'joblib==0.11' 'natsort==5.2.0' 'plotly==2.4.1' 'beakerx==1.2.0' 'pysolr==3.7.0' 'rpy2==2.8.5' \
    'scanpy==0.4.2' 'tqdm==4.23.4' 'validators==0.12.0' 'pysam==0.15.1' 'hca==4.4.10' 'simplegeneric==0.8.1'"

# Create the Python 3.7 environment
RUN conda create -y --name python3.7 python=3.7 anaconda ipykernel==5.1.0

# Install the GenePattern Notebook Environment packages. Running these commands in a new shell isolates
# the conda environment.
RUN /bin/bash -c "source activate python3.7 && \
    conda install -y -c conda-forge 'biopython' 'igraph' 'babel' 'louvain' 'tzlocal' "
RUN /bin/bash -c "source activate python3.7 && \
    conda install -y -c bioconda 'pysam'"
RUN /bin/bash -c "source activate python3.7 && \
    pip install nbtools genepattern-python genepattern-notebook jupyter_wysiwyg \
    'igv-jupyter==0.9.7' 'cyjupyter==0.2.0' 'ccalnoir==2.7' 'cuzcatlan==0.9.3' 'ndex2==1.2.0.*' 'py2cytoscape==0.7.0' \
    'plotly==2.4.1' 'rpy2==2.8.5' 'opencv-python==4.0.0.21' 'hca==4.8.0' 'humanfriendly==4.12.1' scanpy memory_profiler"
RUN /bin/bash -c "source activate python3.7 && \
    pip install -e git://github.com/LankyCyril/genefab.git#egg=genefab \
    --src=/opt/conda/envs/python3.7/lib/python3.7/site-packages/"
