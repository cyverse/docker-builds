FROM r-base:latest
MAINTAINER Upendra Devisetty <upendra@cyverse.org>
LABEL Description "This Dockerfile is for edgeR multifactorial"

# Run updates
RUN apt-get update
RUN apt-get upgrade -y
RUN apt-get install -y git
RUN apt-get install r-base-dev -y
RUN apt-get install libxml2 -y
RUN apt-get install libxml2-dev -y
RUN apt-get -y install libcurl4-gnutls-dev
RUN apt-get install -y libssl-dev

RUN Rscript -e 'source("http://bioconductor.org/biocLite.R"); biocLite("edgeR");'
RUN Rscript -e 'source("http://bioconductor.org/biocLite.R"); biocLite("DESeq2");'
RUN Rscript -e 'source("http://bioconductor.org/biocLite.R"); biocLite("genefilter");'
RUN Rscript -e 'install.packages("devtools", dependencies = TRUE);'
RUN Rscript -e 'devtools::install_github("PF2-pasteur-fr/SARTools")'
RUN Rscript -e 'install.packages("getopt", dependencies = TRUE);'
RUN Rscript -e 'install.packages("knitr", dependencies = TRUE);'
RUN Rscript -e 'install.packages("RColorBrewer", dependencies = TRUE);'

# Add multiple custom functions for Pie_compare, Pie_plot and Bar_compare plots
ENV EDGER1234 https://github.com/upendrak/edgeR_multi.git
RUN git clone $EDGER1234

WORKDIR /edgeR_multi

# change permissions to the wrapper script
ENV BINPATH /usr/bin
RUN chmod +x run_edgeR_multi.r && cp run_edgeR_multi.r $BINPATH
RUN rm -r test_data Dockerfile
RUN mv *.* /

ENTRYPOINT ["run_edgeR_multi.r"]
CMD ["-h"]

# Building and testing
# sudo docker build -t"=rbase/edger-multi:1.0" .
# sudo docker run rbase/edger-multi:1.0 -h
# sudo docker run --rm -v $(pwd):/working-dir -w /working-dir rbase/edger-multi:1.0 --project "test_edgeR" --author "Upendra" --Dir raw --OutDir ./ --target target.txt --features "alignment_not_unique","ambiguous","no_feature","not_aligned","too_low_aQual" --replicates 2 --varInt "group" --condRef "WT" --cpmCutoff 1 --genesele "pairwise" --norm "TMM" --alpha 0.05 --pAdjust "BH" --colors "dodgerblue","orange"
