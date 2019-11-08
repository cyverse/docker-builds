FROM ubuntu:18.04
RUN apt-get update && apt-get install -y python3 python3-pip prodigal hmmer
RUN pip3 install biopython pandas matplotlib seaborn numpy pickle-mixin
RUN pip3 install -U scikit-learn
RUN apt-get -y install git
COPY VIBRANT_run.py /
ENTRYPOINT ["python3", "/VIBRANT_run.py"]

