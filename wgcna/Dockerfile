FROM cyversevice/rstudio-base:3.6.0

RUN Rscript -e 'requireNamespace("BiocManager");' 
RUN Rscript -e 'install.packages("BiocManager");'

RUN R -e 'BiocManager::install(c( "Biobase", "AnnotationDbi", "impute", "GO.db", "org.Mm.eg.db", "preprocessCore"));'

RUN rm -rf /tmp/downloaded_packages/ /tmp/*.rds
RUN Rscript -e 'install.packages("WGCNA");' 

