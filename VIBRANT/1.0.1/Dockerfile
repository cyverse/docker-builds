FROM ubuntu:18.04
RUN apt-get update && apt-get install -y python3 python3-pip prodigal hmmer
RUN pip3 install biopython pandas matplotlib seaborn numpy
RUN pip3 install -U scikit-learn==0.21.3
RUN apt-get -y install git
COPY VIBRANT_run.py /
ENTRYPOINT ["python3", "/VIBRANT_run.py"]

