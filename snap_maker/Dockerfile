FROM quay.io/biocontainers/maker:2.31.10--pl526_16
COPY snap_training.sh /
RUN chmod +x /snap_training.sh
RUN cp /snap_training.sh /usr/bin

ENTRYPOINT [ "snap_training.sh" ]
