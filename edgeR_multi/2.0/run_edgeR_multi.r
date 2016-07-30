#!/usr/bin/Rscript

################################################################################
### R script to compare several conditions with the SARTools and edgeR packages
### Upendra Devisetty
### June 16th, 2016
################################################################################

# Load libraries
library(edgeR)
library(genefilter)
library(devtools)
library(SARTools)
library(getopt)
library(knitr)
library(RColorBrewer)
library(genefilter)


args<-commandArgs(TRUE)

options<-matrix(c('project',  'pn', 1,  "character",    # project name
                  'author', 'au', 1,  "character",      # author of the statistical analysis/report 
                  'Dir',  'r',  2,  "character",          # path to the directory containing raw counts files
                  'rawCounts', 'rc', 2, "character",     # path to the combined raw count file
                  'OutDir',  'w',  1,  "character",     # path to the output file 
                  'target', 't',  1,  "character",      # path to the design/target file
                  'features', 'fe', 2,  "character",    # names of the features to be removed (specific HTSeq-count information and rRNA for example)
                  'varInt', 'v',  2,  "character",      # factor of interest
                  'condRef',  'c',  2,  "character",    # reference biological condition
                  'batch',  'b',  2,  "character",      # blocking factor: NULL (default) or "batch" for example
                  'replicates', 'mi',  2,  "integer",   # minReplicates (smallest number of replicates)
                  'cpmCutoff',  'cc', 2,  "integer",    # counts-per-million cut-off to filter low counts
                  'alpha',  'a',  2,  "double",         # threshold of statistical significance
                  'pAdjust',  'p',  2,  "character",    # p-value adjustment method: "BH" (default) or "BY"
                  'genesele',  'gs', 2,  "character",   # selection of the features in MDSPlot
                  'norm',  'n',  2,  "character",       # normalization method: "TMM" (default), "RLE" (DESeq) or "upperquartile"
                  'colors', 'co', 2,  "character",      # vector of colors of each biological condition on the plots
                  'help',   'h',    0,      "logical"),
                            ncol=4,byrow=TRUE)

ret.opts<-getopt(options,args)

if ( !is.null(ret.opts$help) ) {
  cat(getopt(options, usage=TRUE));
  q(status=1);
}

projectName <- ret.opts$project
author  <-  ret.opts$author
targetFile <- ret.opts$target
rawDir <- ret.opts$Dir
rawCounts <- ret.opts$rawCounts
OutDir <- ret.opts$OutDir
featuresToRemove <- ret.opts$features
varInt  <- ret.opts$varInt
condRef <- ret.opts$condRef
batch <- ret.opts$batch
minReplicates <- ret.opts$replicates
gene.selection <- ret.opts$genesele
cpmCutoff <- ret.opts$cpmCutoff
normalizationMethod <- ret.opts$norm
alpha <- ret.opts$alpha
pAdjustMethod <- ret.opts$pAdjust
col <- ret.opts$colors

# loading target file
target <- loadTargetFile(targetFile=targetFile, varInt=varInt, condRef=condRef, batch=batch)

# Raw counts directory
if (!is.null(ret.opts$Dir)) {
source("/loadCountData.R")
counts <- loadCountData(target=target, rawDir=rawDir, header=TRUE, skip=0, featuresToRemove=featuresToRemove)
}

# Raw counts file
if (!is.null(ret.opts$rawCounts)) {
source("/loadCountDatarc.R")
counts <- loadCountDatarc(target=target, rawCounts=rawCounts, header=TRUE, skip=0, featuresToRemove=featuresToRemove)
}

# # description plots
source("/descriptionPlots.r")
majSequences <- descriptionPlots(counts, n=3, group=target[,varInt], output.file=output.file, col)

# edgeR analysis
source("/run.edgeR.r")
out.edgeR <- run.edgeR(counts, target, varInt, condRef, batch, cpmCutoff, minReplicates, normalizationMethod, pAdjustMethod)

# MDS + clustering
source("/clusterPlot.R")
source("/MDSPlot.R")
source("/heatmap.R")
clusterPlot(group=target[,varInt], output.file="cluster.png")  
MDSPlot(group=target[,varInt], gene.selection, col, output.file="MDS.png")
Heatmap(output.file="heatmap.png")

# summary of the analysis (boxplots, dispersions, export table, nDiffTotal, histograms, MA plot)
source("/summarizeResults.edgeR.r")
summaryResults <- summarizeResults.edgeR(out.edgeR, group=target[,varInt], counts, alpha, col)

# generating HTML report
source("/writeReport.edgeR.r")
writeReport.edgeR(target=target, counts=counts, out.edgeR=out.edgeR, summaryResults=summaryResults,
                  majSequences=majSequences, OutDir=OutDir, projectName=projectName, author=author,
                  targetFile=targetFile, rawDir=rawDir, featuresToRemove=featuresToRemove, varInt=varInt,
                  condRef=condRef, batch=batch, alpha=alpha, pAdjustMethod=pAdjustMethod, colors=colors,
                  gene.selection=gene.selection, normalizationMethod=normalizationMethod)

# End
