FROM genomicsiter/nanodj:latest

COPY jupyter_notebook_config.json /opt/conda/etc/jupyter/jupyter_notebook_config.json

USER root

# Install the icommands, curl, and wget
RUN apt-get update \
    && apt-get install -y lsb wget gnupg apt-transport-https python3 python-requests curl \
    && apt-get clean \
    && rm -rf /usr/lib/apt/lists/* \
    && fix-permissions $CONDA_DIR

RUN wget -qO - https://packages.irods.org/irods-signing-key.asc | apt-key add - \
    && echo "deb [arch=amd64] https://packages.irods.org/apt/ xenial main" > /etc/apt/sources.list.d/renci-irods.list \
    && apt-get update \
    && apt-get install -y irods-icommands \
    && apt-get clean \
    && rm -rf /usr/lib/apt/lists/* \
    && fix-permissions $CONDA_DIR


USER jovyan

# install foundational jupyter lab
RUN conda update -n base conda \
    && conda install jupyterlab=0.35.4 \
    && conda clean -tipsy \
    && fix-permissions $CONDA_DIR

# install jupyter hub and extra doodads
RUN jupyter lab --version \
    && jupyter labextension install @jupyterlab/hub-extension@0.12.0 \
                                    @jupyter-widgets/jupyterlab-manager@0.38.1 \
                                    jupyterlab_bokeh@0.6.3

# Install the irods plugin for jupyter lab
RUN pip install jupyterlab_irods==2.0.1 \
    && jupyter serverextension enable --py jupyterlab_irods \
    && jupyter labextension install ijab@2.0.1

COPY entry.sh /bin
RUN mkdir /home/jovyan/.irods

ENTRYPOINT ["bash", "/bin/entry.sh"]
