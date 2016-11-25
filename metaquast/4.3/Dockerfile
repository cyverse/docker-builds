FROM ubuntu:14.04.3
MAINTAINER Upendra Devisetty <upendra@cyverse.org>
LABEL Description "This Dockerfile is for metaQUAST-4.3"

RUN apt-get update 

# Dependencies
RUN apt-get install -y g++ \
                          make \
                          wget \
                          python \
                          python-matplotlib \
                          zlib1g-dev \
                          cmake \
                          openjdk-6-jdk \
                          curl \
                          libboost-all-dev \   
                          libncurses5-dev

# cpanm modules
RUN curl -L http://cpanmin.us | perl - App::cpanminus
RUN cpanm Time::HiRes

# Download Metaquast-4.3
RUN wget https://downloads.sourceforge.net/project/quast/quast-4.3.tar.gz
RUN tar -xzf quast-4.3.tar.gz

# Specify entrypoint
ENTRYPOINT ["/quast-4.3/metaquast.py"]
CMD ["-h"]

# Build the image
# docker build -t"=ubuntu/metaquast:4.3" .
# Testing the image
# Without any arguments
# sudo docker run ubuntu/metaquast-4.3 -h
# With arguments
# Test dataset
# sudo docker run ubuntu/metaquast-4.3 --test
# Metaquast with ref
# sudo docker run --rm -v $(pwd):/working-dir -w /working-dir ubuntu/metaquast-4.0 -o metaquast_test_output -R meta_ref_1.fasta,meta_ref_2.fasta,meta_ref_3.fasta meta_contigs_1.fasta meta_contigs_2.fasta
# Metaquast with no reference
# sudo docker run --rm -v $(pwd):/working-dir -w /working-dir ubuntu/metaquast-4.0 -o metaquast_test_output_no_ref meta_contigs_1.fasta meta_contigs_2.fasta



