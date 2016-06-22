FROM r-base:latest
MAINTAINER Upendra Devisetty <upendra@cyverse.org>
LABEL Description "This Dockerfile is for IUTA"

# Run updates
RUN apt-get update
RUN apt-get upgrade -y
RUN apt-get install wget -y
RUN apt-get install r-base-dev -y
RUN apt-get install libxml2 -y
RUN apt-get install libxml2-dev -y

RUN wget http://www.niehs.nih.gov/resources/files/IUTA_1.0.tar.gz
RUN tar zxvf IUTA_1.0.tar.gz

RUN Rscript -e 'source("http://bioconductor.org/biocLite.R"); biocLite("Rsamtools");'
RUN Rscript -e 'install.packages("/IUTA_1.0.tar.gz", repos = NULL, type="source");'
RUN Rscript -e 'install.packages("getopt");'
RUN Rscript -e 'install.packages("ggplot2");'
RUN Rscript -e 'install.packages("reshape2");'
RUN Rscript -e 'install.packages("grid");'

# Add multiple custom functions for Pie_compare, Pie_plot and Bar_compare plots
ADD pie_compare.R /
ADD pie_plot.R /
ADD bar_compare.R /

# Add wrapper script
ADD run_IUTA.R /
RUN chmod +x /run_IUTA.R && cp /run_IUTA.R /usr/bin

ENTRYPOINT ["run_IUTA.R"]
CMD ["-h"]

# Building and testing
# sudo docker build -t"=ubuntu/iuta:1.0" .
# sudo docker run ubuntu/iuta:1.0 -h
# sudo docker run --rm -v $(pwd):/working-dir -w /working-dir ubuntu/iuta:1.0 --gtf mm10_kg_sample_IUTA.gtf --bam1 bam_1 --bam2 bam_2 --fld empirical --test.type SKK,CQ,KY --numsamp 3 --output IUTA_test_1 --groups 1,2 --gene.id Pcmtd1 
# With no gene id (all genes compressed)
# sudo docker run --rm -v $(pwd):/working-dir -w /working-dir ubuntu/iuta:1.0 --gtf mm10_kg_sample_IUTA.gtf --bam1 bam_1 --bam2 bam_2 --fld empirical --test.type SKK,CQ,KY --numsamp 3 --output IUTA_test_1 --groups 1,2
