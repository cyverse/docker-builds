FROM cyversevice/jupyterlab-scipy:latest

USER root

RUN pip install genepattern-notebook && \
    jupyter nbextension enable --py widgetsnbextension && \
    jupyter nbextension install --py genepattern && \
    jupyter nbextension enable --py genepattern && \
    jupyter serverextension enable --py genepattern && \
    jupyter nbextension install --py nbtools && \
    jupyter nbextension enable --py nbtools && \
    jupyter serverextension enable --py nbtools

USER jovyan
