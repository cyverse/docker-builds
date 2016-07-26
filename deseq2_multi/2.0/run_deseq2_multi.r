#!/usr/bin/Rscript

################################################################################
### R script to compare several conditions with the SARTools and DESeq2 packages
### Upendra Devisetty
### June 16th, 2016
################################################################################

# Load libraries
library(DESeq2)
library(genefilter)
library(devtools)
library(SARTools)
library(getopt)
library(knitr)
library(RColorBrewer)
library(genefilter)


args<-commandArgs(TRUE)

options<-matrix(c('project',  'pn', 1,  "character",      # project name
                  'author', 'au', 1,  "character",        # author of the statistical analysis/report  
                  'Dir',  'r',  2,  "character",          # path to the directory containing raw counts files
                  'rawCounts', 'rc', 2, "character",      # path to the combined raw count file
                  'OutDir',  'w',  1,  "character",       # path to the output file 
                  'target', 't',  1,  "character",        # path to the design/target file
                  'features', 'fe', 2,  "character",      # names of the features to be removed (specific HTSeq-count information and rRNA for example)
                  'varInt', 'v',  2,  "character",        # factor of interest
                  'condRef',  'c',  2,  "character",      # reference biological condition
                  'batch',  'b',  2,  "character",        # blocking factor: NULL (default) or "batch" for example
                  'locfunc', 'l', 2, "character",         # "median" (default) or "shorth" to estimate the size factors
                  'fitType', 'f', 2, "character",         # mean-variance relationship: "parametric" (default) or "local"
                  'cooksCutoff', 'cc',  2, "logical",     # outliers detection threshold (TRUE to let DESeq2 choosing it or FALSE to disable the outliers detection)
                  'independentFiltering', 'if', 2, "logical",   # TRUE or FALSE to perform the independent filtering or not
                  'typeTrans', 'tt', 2, "character",    # transformation method for PCA/clustering with DESeq2: "VST" or "rlog"
                  'alpha',  'a',  2,  "double",           # threshold of statistical significance
                  'pAdjust',  'p',  2,  "character",      # p-value adjustment method: "BH" (default) or "BY"
                  'colors', 'co', 2,  "character",        # vector of colors of each biological condition on the plots
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
fitType <- ret.opts$filtType
cooksCutoff <- ret.opts$cooksCutoff
independentFiltering <- ret.opts$independentFiltering
locfunc <- ret.opts$locfunc
alpha <- ret.opts$alpha
pAdjustMethod <- ret.opts$pAdjust
typeTrans <- ret.opts$typeTrans
col <- ret.opts$colors


# loading target file
target <- loadTargetFile(targetFile=targetFile, varInt=varInt, condRef=condRef, batch=batch)

# loading counts

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

# description Plots
source("/descriptionPlots.r")
majSequences <- descriptionPlots(counts=counts, n=3, group=target[,varInt], output.file=output.file, col=col)

# analysis with DESeq2
source("/run.DESeq2.r")
out.DESeq2 <- run.DESeq2(counts=counts, target=target, varInt=varInt, batch=batch, locfunc=locfunc, fitType=fitType, pAdjustMethod=pAdjustMethod,
                         cooksCutoff=cooksCutoff, independentFiltering=independentFiltering, alpha=alpha)
# MDS + clustering
source("/exploreCounts.R")
exploreCounts(object=out.DESeq2$dds, group=target[,varInt], typeTrans=typeTrans, col=col)

# summary of the analysis (boxplots, dispersions, export table, nDiffTotal, histograms, MA plot)
source("/summarizeResults.DESeq2.r")
summaryResults <- summarizeResults.DESeq2(out.DESeq2, group=target[,varInt], col=col, independentFiltering=independentFiltering,
                                          cooksCutoff=cooksCutoff, alpha=alpha)

# generating HTML report
source("/writeReport.DESeq2.r")
writeReport.DESeq2(target=target, counts=counts, out.DESeq2=out.DESeq2, summaryResults=summaryResults, majSequences=majSequences, OutDir=OutDir, 
		   projectName=projectName, author=author, targetFile=targetFile, rawDir=rawDir, featuresToRemove=featuresToRemove, varInt=varInt,
       condRef=condRef, batch=batch, fitType=fitType, cooksCutoff=cooksCutoff, independentFiltering=independentFiltering, alpha=alpha, 
		   pAdjustMethod=pAdjustMethod, typeTrans=typeTrans, locfunc=locfunc, colors=col)
# End
