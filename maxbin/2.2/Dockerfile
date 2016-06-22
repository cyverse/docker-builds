FROM ubuntu:14.04.3
MAINTAINER Upendra Kumar Devisetty <upendra@cyverse.org>
LABEL Description="This Dockerfile is used for building Maxbin-2.2 Docker image" 
# To get rid of all the messages
RUN DEBIAN_FRONTEND=noninteractive
# To update the image
RUN apt-get update
# To install all the dependencies
RUN apt-get install -y build-essential wget make curl unzip python
# To download the Maxbin software and untar it
RUN wget https://sourceforge.net/projects/maxbin/files/latest/download
RUN tar xvf download
# To se the Workdirectory
WORKDIR /MaxBin-2.2
# To make
RUN cd /MaxBin-2.2/src && make
# To install the dependencies for Maxbin
RUN /MaxBin-2.2/autobuild_auxiliary
# To change permission of the Maxbin script
RUN chmod +x /MaxBin-2.2/run_MaxBin.pl
# Entrypoint
ENTRYPOINT ["/MaxBin-2.2/run_MaxBin.pl"]

# TO build the image
# sudo docker build -t="ubuntu/maxbin:2.2" .

# To test the image without arguments
# sudo docker run ubuntu/maxbin:2.2

# To test the image with arguments
# mkdir -p docker_test
# cd docker_test
# Download the test files
# wget http://downloads.jbei.org/data/microbial_communities/MaxBin/getfile.php?20x.scaffold
# mv getfile.php\?20x.scaffold 20x.scaffold
# wget http://downloads.jbei.org/data/microbial_communities/MaxBin/getfile.php?20x.abund
# mv getfile.php\?20x.abund 20x.abund
# sudo docker run -v /home/upendra_35/MaxBin-2.2/docker_test:/binning -w /binning ubuntu/maxbin:2.2 -contig 20x.scaffold -abund 20x.abund -out 20x.out -thread 4

