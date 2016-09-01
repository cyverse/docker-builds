FROM ubuntu:14.04.3
MAINTAINER Anne Bowen <adb@tacc.utexas.edu>
LABEL Description "This Dockerfile is for creating heatmap"

RUN apt-get update
RUN apt-get -y  install  build-essential python python-dev python-pip python-setuptools python-tk gfortran libblas-dev libxft-dev libfreetype6 libfreetype6-dev libpng-dev libatlas-dev liblapack-dev
RUN easy_install pip
RUN pip install --upgrade distribute
RUN pip install numpy==1.7.1
RUN pip install scipy==0.12.0
RUN pip install matplotlib
RUN pip install pandas
RUN pip install names
COPY heatmap.py /usr/bin
ENTRYPOINT ["python","/usr/bin/heatmap.py"]

# Build
#sudo docker build -t"=ubuntu/heatmap:1.0" .

# Running without any arguments
#sudo docker run ubuntu/heatmap:1.0

# Testing
#sudo docker run --rm -v ~/heatmap/test_data:/heatmap -w /heatmap ubuntu/heatmap:1.0 --i data.txt
#sudo docker run --rm -v ~/heatmap/test_data:/heatmap -w /heatmap ubuntu/heatmap:1.0 --i data.txt --row_header "complete" --column_method "weighted"