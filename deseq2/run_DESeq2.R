#!/usr/bin/env Rscript

## Example command line invocation:
## "./Rscript.exe" "./run_DESeq2.R" --input "./DESeq_test_data.tsv" --namecol 1 --condition "untreated,untreated,untreated,untreated,treated,treated,treated" --libType "single-end,single-end,paired-end,paired-end,single-end,paired-end,paired-end" --pair "untreated,treated" --fdr 0.1 --filter 0.4

## Example input file (pasilla data set)
## 
## gene_id	untreated1	untreated2	untreated3	untreated4	treated1	treated2	treated3
## FBgn0000003	0	0	0	0	0	0	1
## FBgn0000008	92	161	76	70	140	88	70
## FBgn0000014	5	1	0	0	4	0	0
## FBgn0000015	0	2	1	2	1	0	0
## FBgn0000017	4664	8714	3564	3150	6205	3072	3334

# Auto-set up dependencies
# This is for API deployment
#
#  .libPaths(getwd())
#  install.packages("getopt", repos="http://lib.stat.cmu.edu/R/CRAN", lib=getwd())
#  install.packages("gplots", repos="http://lib.stat.cmu.edu/R/CRAN", lib=getwd())
#  source("http://bioconductor.org/biocLite.R")
#  biocLite("DESeq2", destdir=getwd(), lib=getwd())
# .libPaths(getwd())

library("getopt")
library("DESeq2")
library("RColorBrewer")
library("gplots")

args<-commandArgs(TRUE)

#############################################################
## Section 0: Command-line preparations and getopt parsing ##
#############################################################

##
## Creates a matrix called options. We're going to need 
## several different options for DESeq deployment: input, 
## namecol, condition, libType, pair, fdr, filter, and 
## mtc (holm, hochberg, hommel, bonferroni, BH, BY, fdr, none).
##
options<-matrix(c('input',		'i',	1,	"character",
			'namecol',		'n',	1,	"integer",
			'condition',	'c',	1,	"character",
			'libType',		'l',	1,	"character",
			'pair',		'p',	1,	"character",
			'fdr',		'f',	1,	"double",
			'filter',		'x',	1,	"double",
			'mtc',		'm',	0,	"character"),
			ncol=4,byrow=TRUE)

ret.opts<-getopt(options,args)

##
## Details how the data file is formatted.
##
tableDelim <- "\t"
tableQuote <- "\""
tableDec <- "."
tableComment <- "#"
tableHeader <- FALSE

##
## == REQUIRED ==
## tableFile is the location of the input file.
## tableRowNames corresponds to the index of
## the column that contains the gene names.
##
tableFile <- ret.opts$input
tableRowNames <- ret.opts$namecol

##
## == REQUIRED ==
## condition is a single, comma-delimited string containing the info 
## about what conditions/factors are present in the data file. For 
## example, the pasilla conditions look like this:
## "untreated,untreated,untreated,untreated,treated,treated,treated"
## 
## libType is a single, comma-delimited string containing all the 
## info about the library types (single-end or paired end) of each 
## condition. For example, the pasilla libTypes looks like this: 
## "single-end,single-end,paired-end,paired-end,single-end,paired-end,paired-end"
##
condition <- unlist(strsplit(ret.opts$condition, ","))
libType <- unlist(strsplit(ret.opts$libType, ","))

##
## == REQUIRED ==
## pairFactor corresponds to the comma-separated
## pair of factors used for comparison.
##
pairFactor <- unlist(strsplit(ret.opts$pair, ","))

##
## == REQUIRED ==
## minFDR corresponds to the minimum false-discovery 
## rate. The DESeq vignette uses 0.1 this value.
##
minFDR <- ret.opts$fdr

##
## == REQUIRED ==
## Filter is used to remove genes with
## little statistical significance.
##
filter <- ret.opts$filter

##
## == RECOMMENDED ==
## mtcMethod corresponds to the user-specified
## method for multiple testing correction (holm, 
## hochberg, hommel, bonferroni, BH, BY, fdr, none).
##
if(is.null(ret.opts$mtc))
{
	mtcMethod <- "BH"
} else
{
	mtcMethod <- ret.opts$mtc
}

############################################
## Section 1: Input data and preparations ##
############################################

##
## Read the data file.
##
originalCountTable <- read.table(tableFile, row.names=tableRowNames, 
	header = TRUE, sep=tableDelim)
head(originalCountTable)

##
## groupFactor corresponds to the names 
## of the columns in the input file.
##
groupFactor <- colnames(originalCountTable)
gf <- factor(unlist(strsplit(groupFactor, ",")))
group <- factor(gf)

##
## Read the metadata file. If a metadata file was not
## provided, attempt to create one from scratch.
##
originalDesign = data.frame(row.names = colnames(originalCountTable), condition, libType)
originalDesign

##
## Determine how many replicates there are.
##
conditionSummary <- table(condition)
replicates <- "all"
none <- TRUE
for(i in 1:length(conditionSummary))
{
	if(conditionSummary[i] == 1)
	{
		if(replicates == "all")
		{
			replicates <- "some"
		}
	} else
	{
		none <- FALSE
	}
}
if(none == TRUE)
{
	replicates <- "none"
}

####################################################
## Section 2: Dispersion Estimation and Filtering ##
####################################################

##
## Instantiate a newCountDataSet, the central
## data structure in the DESeq package.
##

cdsFull<-DESeqDataSetFromMatrix(originalCountTable,DataFrame(condition),~condition)
#cdsFull = newCountDataSet(originalCountTable, originalDesign)

if(replicates != "none"){
	##
	## Normalisation step (no relation to normality/normal 
	## distribution). Brings count values to a common scale.
	##
	cdsFull = estimateSizeFactors(cdsFull)
	sizeFactors(cdsFull)
	head(counts(cdsFull, normalized=TRUE))
	##
	## Estimate dispersions and plot them.
	##
	cdsFull = estimateDispersions(cdsFull)
	png(paste("DESeq_Dispersions",".png",sep=""))
	plotDispEsts(cdsFull)
	dev.off()
}else{
	##
	## Estimate dispersions blindly and plot them.
	##
	cdsFull = estimateDispersions(cdsFull, method="blind", sharingMode="fit-only")
	png(paste("DESeq_Dispersions",".png",sep=""))
	plotDispEsts(cdsFull)
	dev.off()
}

##
## Consider a filter criterion rs, the overall sum of 
## counts (irrespective of biological condition), and
## remove the genes in a user-supplied lower quantile.
##
rs = rowSums(counts(cdsFull))
use = (rs > quantile(rs, probs=filter))
table(use)
cdsFilt = cdsFull[use,]
stopifnot(!any(is.na(use)))


###########################################################
## Section 3: Inference: Calling differential expression ##
###########################################################

cdsFull<-nbinomWaldTest(cdsFull);
dds<-DESeq(cdsFull)
resFull<-results(dds)
head(resFull)

## Plot the log2 fold changes against the mean
## normalised counts, coloring in red those
## genes that are signicant at 10% FDR, and
## output a histogram of the p-values.
##
png(paste("DESeq_MAplot",".png",sep=""))
plotMA(resFull)
dev.off()
png(paste("DESeq_pValues",".png",sep=""))
hist(resFull$pval, breaks=100, col="skyblue", border="slateblue", main="")
dev.off()

##
## Filter for significant genes, according to
## a user-supplied False Discovery Rate.
##
resFullSig = resFull[!is.na(resFull$padj) < minFDR,]

##
## List the most signicantly differentially expressed
## genes, as well as the most strongly down-regulated
## and up-regulated of the signicant genes.
##
#head(resFullSig[order(resFullSig$pval),])
#head(resFullSig[order(resFullSig$foldChange, -resFullSig$baseMean),])
#head(resFullSig[order(-resFullSig$foldChange, -resFullSig$baseMean),])

##
## Write results to tab-delimited files.
##
write.table(resFull, file=paste("DESeq_results", "txt", sep="."),
sep="\t", quote=FALSE)
write.table(resFullSig, file=paste("DESeq_results_significant", "txt", sep="."),
sep="\t", quote=FALSE)

#Now we want to transform the raw discretely
#distributed counts so that we can do clustering.

rld <- rlogTransformation(dds, blind=TRUE)
vsd <- varianceStabilizingTransformation(dds, blind=TRUE)

# good heatmap!

select <- order(rowMeans(counts(dds,normalized=TRUE)),decreasing=TRUE)[1:30]
hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)
heatmap.2(counts(dds,normalized=TRUE)[select,], col = hmcol,
Rowv = FALSE, Colv = FALSE, scale="none",
dendrogram="none", trace="none", margin=c(10,6))
#dev.copy(png,"DESeq2_heatmap1")
png(paste("DESeq2_heatmap1",".png",sep=""))
dev.off()
heatmap.2(assay(rld)[select,], col = hmcol,
Rowv = FALSE, Colv = FALSE, scale="none",
dendrogram="none", trace="none", margin=c(10, 6))
#dev.copy(png,"DESeq2_heatmap2")
png(paste("DESeq2_heatmap2",".png",sep=""))
dev.off()
heatmap.2(assay(vsd)[select,], col = hmcol,
Rowv = FALSE, Colv = FALSE, scale="none",
dendrogram="none", trace="none", margin=c(10, 6))
#dev.copy(png,"DESeq2_heatmap3")
png(paste("DESeq2_heatmap3",".png",sep=""))
dev.off()

##
## List of files that have been created:
##
## DESeq_Dispersion.png
## DESeq_MAplot.png
## DESeq_pValues.png
## DESeq_results.txt
## DESeq_results_significant.txt
##
