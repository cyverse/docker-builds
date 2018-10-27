#!/usr/local/bin/Rscript

# Load libraries
library(oposSOM)
library(biomaRt)
library(stringi)
library(getopt)

args<-commandArgs(TRUE)

options<-matrix(c('file',  'f', 1,  "character",       
                  'datasetname', 'n', 1,  "character",       
                  'databasebiomart',  'b',  1,  "character",
                  'databasehost',  'ho',  1,  "character",       
                  'databasedataset', 'd', 1, "character", 
                  'databaseidtype', 'i', 1, "character",  
                  'samples',  'l',  1,  "character",
                  'log10', 'g', 0, "logical",    
                  'help',   'h',    0,      "logical"),
                   ncol=4,byrow=TRUE)

ret.opts<-getopt(options,args)

if ( !is.null(ret.opts$help) ) {
  cat(getopt(options, usage=TRUE));
  q(status=1);
}

file <- ret.opts$file
datasetname  <-  ret.opts$datasetname
databasebiomart <- ret.opts$databasebiomart
databasehost <- ret.opts$databasehost
databasedataset <- ret.opts$databasedataset
databaseidtype <- ret.opts$databaseidtype
file2 <- ret.opts$samples

data <- read.csv(file = file, sep = "\t")

env <- opossom.new()

env$preferences <- list(dataset.name = datasetname, dim.1stLvlSom = "auto",
                       dim.2ndLvlSom = 20, training.extension = 1, rotate.SOM.portraits = 0,
                       flip.SOM.portraits = FALSE, activated.modules = list(reporting = TRUE,
                                                                            primary.analysis = TRUE, sample.similarity.analysis = TRUE,
                                                                            geneset.analysis = TRUE, geneset.analysis.exact = TRUE,
                                                                            group.analysis = TRUE, difference.analysis = TRUE),
                       database.biomart = databasebiomart, database.host = databasehost,
                       database.dataset = databasedataset, database.id.type = databaseidtype, standard.spot.modules = "dmap",
                       spot.coresize.modules = 3, spot.threshold.modules = 0.95,
                       spot.coresize.groupmap = 5, spot.threshold.groupmap = 0.75,
                       adjust.autogroup.number = 0, feature.centralization = TRUE,
                       sample.quantile.normalization = TRUE, pairwise.comparison.list = NULL)

# If log transformation 

if ( !is.null(ret.opts$log10) ) {
  data_log <- log10(data[,-1] + 1);
  data_log <- cbind(data[,1],data_log);
  names(data_log)[1] <- names(data[1]);
  print("log transformation");
  env$indata <- data_log
} else  {
# If not log transformation
  print("Not log transformation");
  env$indata <- data
}

# Sample file (txt file containing sample lables and colors)

oposom <- read.csv(file2, sep = "\t")
oposom$label <- as.character(oposom$label)
oposom$color <- as.character(oposom$color)

vector <- c()
for (i in 1:nrow(oposom)){
  vector <- c(vector, rep(oposom$label[i], oposom$replicates[i]))
}            
env$group.labels <- vector
          
vector2 <- c()
for (i in 1:nrow(oposom)){
  vector2 <- c(vector2, rep(oposom$color[i], oposom$replicates[i]))
}            
env$group.colors <- vector2

# Run opossom

opossom.run(env)
