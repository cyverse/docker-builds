FROM python:2.7.10-slim

WORKDIR /root

COPY ht_path_list_splitter.py ./
COPY config config/

ENTRYPOINT [ "python", "/root/ht_path_list_splitter.py" ]
CMD [ "--help" ]
