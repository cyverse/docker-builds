FROM cyverse/jupyterlab-qiime2:2018.11

RUN pip install bash_kernel
RUN python -m bash_kernel.install

COPY entry.sh /bin
RUN mkdir /home/jovyan/.irods

ENTRYPOINT ["bash", "/bin/entry.sh"]

