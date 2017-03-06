## Dockerfile
FROM ubuntu:14.04.3
MAINTAINER Amanda Cooksey	
LABEL Description="cufflinks-2.2.1.Linux update"

# Install all the updates and download dependencies
RUN apt-get update && apt-get install -y git wget

# Download the cufflinks-2.2.1
RUN wget -O- http://cole-trapnell-lab.github.io/cufflinks/assets/downloads/cufflinks-2.2.1.Linux_x86_64.tar.gz | tar xzvf -
ENV PATH /cufflinks-2.2.1.Linux_x86_64/:$PATH

# ADD the cuffcompare-2.2.1 wrapper script
ADD cuffcompare-2.2.1.pl /usr/bin

# Change the permissions and the path for the wrapper script
RUN chmod +x /usr/bin/cuffcompare-2.2.1.pl

# Entrypoint
ENTRYPOINT ["/usr/bin/cuffcompare-2.2.1.pl"]

# Build dockerimage from dockerfile
# sudo docker build -t"=ubuntu/cuffcompare-2.2.1:latest" .

# cd to sample_files directory and then do test run
# sudo docker run --rm -v $(pwd):/home/amcooksey/cuffcompare-2.2.1/sample_files -w /home/amcooksey/cuffcompare-2.2.1/sample_files ubuntu/cuffcompare-2.2.1:latest --infile flower4_transcripts.gtf --infile flower6-7_transcripts.gtf 

# Push the dockerfile to the gitrepo

