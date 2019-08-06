FROM quay.io/biocontainers/parallel-fastq-dump:0.6.5--py_0
MAINTAINER Upendra Devisetty <upendra@cyverse.org>
LABEL Description "This Dockerfile is for parallel-fastq-dump multiple files"

ADD parallel_fastq_dump_wrapper.sh /
RUN chmod +x /parallel_fastq_dump_wrapper.sh && mv /parallel_fastq_dump_wrapper.sh /usr/bin

ENTRYPOINT ["parallel_fastq_dump_wrapper.sh"] 
CMD ["-h"]
