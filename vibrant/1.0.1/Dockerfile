from ubuntu:18.04
RUN apt-get update && apt-get install -y python3 python3-pip prodigal hmmer
RUN pip3 install biopython pandas matplotlib seaborn numpy pickle-mixin
RUN pip3 install -U scikit-learn
RUN apt-get -y install git
RUN git clone https://github.com/AnantharamanLab/VIBRANT
RUN chmod -R 777 VIBRANT
WORKDIR VIBRANT
RUN cp VIBRANT_run.py /usr/bin

ENTRYPOINT ["python3","/usr/bin/VIBRANT_run.py"]
