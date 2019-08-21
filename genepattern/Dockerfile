FROM genepattern/genepattern-notebook:19.08

USER root

RUN apt-get update \
        && apt-get install -y jq \
        && cd /opt/conda/etc/jupyter \
        && mv jupyter_notebook_config.json jupyter_notebook_config.json.bak \
        && jq '.NotebookApp = .NotebookApp + {"token":"","password":""}' jupyter_notebook_config.json.bak \
            > jupyter_notebook_config.json \
        && rm jupyter_notebook_config.json.bak

USER $NB_USER
