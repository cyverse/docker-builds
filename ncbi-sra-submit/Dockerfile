FROM python:2.7-slim

WORKDIR /root

COPY requirements.txt ./
RUN set -x \
    && apt-get update \
    && apt-get install -y gcc libxml2-dev libxslt1-dev lib32z1-dev --no-install-recommends \
    && rm -rf /var/lib/apt/lists/* \
    && pip install -r requirements.txt \
    && apt-get purge -y --auto-remove gcc lib32z1-dev

# Download Aspera Connect client from http://downloads.asperasoft.com/connect2/
ADD http://download.asperasoft.com/download/sw/connect/3.5/aspera-connect-3.5.1.92523-linux-64.sh aspera-connect-install.sh

# Install Aspera Connect client to ~/.aspera
RUN chmod 755 aspera-connect-install.sh
RUN ./aspera-connect-install.sh \
    && rm aspera-connect-install.sh

COPY ncbi_sra_submit.py metadata_client.py ncbi_sra_report_download.py ./

VOLUME [ "/root/config", "/root/templates", "/root/schemas" ]

ENTRYPOINT [ "python", "/root/ncbi_sra_submit.py" ]
CMD [ "--help" ]

# Example run commands:
# docker run --rm -e IPLANT_USER -e IPLANT_EXECUTION_ID -w /de-app-work -v /condor/scratch/workdir:/de-app-work --volumes-from=ncbi-ssh-key:ro --volumes-from=ncbi-sra-configs:ro --net=bridge discoenv/ncbi-sra-submit --submit-mode create --input-metadata metadata.json --input-dir BioProject
# docker run --rm -e IPLANT_USER -e IPLANT_EXECUTION_ID -w /de-app-work -v /condor/scratch/workdir:/de-app-work --volumes-from=ncbi-ssh-key:ro --volumes-from=ncbi-sra-configs:ro --net=bridge --entrypoint /root/ncbi_sra_report_download.py discoenv/ncbi-sra-submit --submit-dir 'username.bio-project-folder-id'
