FROM python:3-slim

WORKDIR /root

RUN set -x \
    && apt-get update \
    && apt-get install -y gcc libxml2-dev libxslt1-dev lib32z1-dev --no-install-recommends \
    && rm -rf /var/lib/apt/lists/* \
    && apt-get purge -y --auto-remove gcc lib32z1-dev

# Aspera Connect does not allow the root user to install its tool.
RUN useradd --create-home --shell /bin/bash cyverse
USER cyverse
WORKDIR /home/cyverse

COPY requirements.txt ./
RUN pip install -r requirements.txt

# Download Aspera Connect client from https://www.ibm.com/aspera/connect/
# and install client to ~/.aspera
ADD --chown=cyverse https://d3gcli72yxqn2z.cloudfront.net/downloads/connect/latest/bin/ibm-aspera-connect_4.2.6.393_linux_x86_64.tar.gz ./
RUN tar -xvf ibm-aspera-connect_4.2.6.393_linux_x86_64.tar.gz \
    && rm ibm-aspera-connect_4.2.6.393_linux_x86_64.tar.gz \
    && mv ibm-aspera-connect_4.2.6.393_linux_x86_64.sh aspera-connect-install.sh \
    && chmod 755 aspera-connect-install.sh \
    && ./aspera-connect-install.sh \
    && rm aspera-connect-install.sh

COPY bioproject.py metadata_client.py ncbi_report_download.py ncbi_submit.py ./

VOLUME [ "/home/cyverse/config", "/home/cyverse/templates", "/home/cyverse/schemas" ]

ENTRYPOINT [ "python", "/home/cyverse/ncbi_submit.py" ]
CMD [ "--help" ]
