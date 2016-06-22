FROM ubuntu:14.04.3
MAINTAINER Upendra Devisetty <upendra@cyverse.org>
LABEL Description "This Dockerfile is for QUAST-4.0"

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

# Clone Quast git repo
RUN wget -O- https://github.com/ablab/quast/archive/release_4.0.tar.gz | tar zxvf -
RUN chmod +x /quast-release_4.0/quast.py

# configure samtools
WORKDIR /quast-release_4.0/libs/samtools
RUN ./configure && make && make install

# Specify entrypoint
ENTRYPOINT ["/quast-release_4.0/quast.py"]
CMD ["-h"]

# Build the image
# docker build -t"=ubuntu/quast-4.0" dockerfile-quast
# Testing the image
# Without any arguments
# sudo docker run ubuntu/quast-4.0 -h
# With arguments
# sudo docker run --rm -v $(pwd):/working-dir -w /working-dir ubuntu/quast-4.0 contigs_1.fasta contigs_2.fasta -R reference.fasta.gz -O operons.txt -G genes.txt -o quast_test_out
# SV calling
# sudo docker run --rm -v $(pwd):/working-dir -w /working-dir ubuntu/quast-4.0 -o quast_test_output_sv -R reference.fasta.gz -O operons.gff -G genes.gff --gage  --gene-finding  --eukaryote  --glimmer  -1 reads1.fastq.gz -2 reads2.fastq.gz contigs_1.fasta contigs_2.fasta

