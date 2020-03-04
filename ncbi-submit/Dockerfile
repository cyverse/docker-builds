FROM python:3-slim

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

COPY bioproject.py metadata_client.py ncbi_report_download.py ncbi_submit.py ./

VOLUME [ "/root/config", "/root/templates", "/root/schemas" ]

ENTRYPOINT [ "python", "/root/ncbi_submit.py" ]
CMD [ "--help" ]
