FROM discoenv/jupyter-lab:beta

USER root

# Install the icommands, curl, and wget
RUN apt-get update \
    && apt-get install -y lsb wget gnupg apt-transport-https python3.6 python-requests curl \
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

# Download qiime2 yml file
RUN wget https://data.qiime2.org/distro/core/qiime2-2018.11-py35-linux-conda.yml

# install foundational jupyter lab
RUN conda env create -n qiime2-2018.11 --file qiime2-2018.11-py35-linux-conda.yml \
    && conda clean -tipsy \
    && rm qiime2-2018.11-py35-linux-conda.yml

WORKDIR /home/jovyan

ENV PATH /opt/conda/envs/qiime2-2018.11/bin/:$PATH

ENTRYPOINT ["jupyter"]

CMD ["lab", "--no-browser"] 
