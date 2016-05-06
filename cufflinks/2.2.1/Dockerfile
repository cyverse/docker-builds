## Dockerfile
FROM ubuntu:14.04.3
MAINTAINER Upendra Devisetty
LABEL Description="cufflinks-2.2.1.Linux update"

# Install all the updates and download dependencies
RUN apt-get update && apt-get install -y git wget

# Download the cufflinks-2.2.1
RUN wget -O- http://cole-trapnell-lab.github.io/cufflinks/assets/downloads/cufflinks-2.2.1.Linux_x86_64.tar.gz | tar xzvf -

# Clone the cufflinks-2.2.1 wrapper script
ENV CUFF2GIT https://github.com/upendrak/cufflinks-2.2.1.git
RUN git clone $CUFF2GIT

# Change the permissions and the path for the wrapper script
RUN chmod +x /cufflinks-2.2.1/cufflinks-2.2.1.pl && cp /cufflinks-2.2.1/cufflinks-2.2.1.pl /usr/bin 
ENV PATH /cufflinks-2.2.1.Linux_x86_64/:$PATH

# Entrypoint
ENTRYPOINT ["/usr/bin/cufflinks-2.2.1.pl"]
CMD ["-h"]

# Build dockerimage from dockerfile
# sudo docker build -t"=ubuntu/cufflinks-2.2.1:latest" .

# cd to sample_files directory and then do test run
# sudo docker run --rm -v $(pwd):/home/upendra_35/cufflinks-2.2.1/sample_files -w /home/upendra_35/cufflinks-2.2.1/sample_files ubuntu/cufflinks-2.2.1:latest --infile SRR070570_WT.fastq.tophat.bam 

# Push the dockerfile to the gitrepo
