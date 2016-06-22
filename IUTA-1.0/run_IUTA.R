#!/usr/bin/Rscript

# Install dependencies
library(Rsamtools)
library(IUTA)
library(getopt)
library(ggplot2)
library(reshape2)
library(grid)


args<-commandArgs(TRUE)

#############################################################
## Command-line preparations and getopt parsing ##
#############################################################

options<-matrix(c('gtf',	'i',	1,	"character",
		  'bam1',	'ba1',	1,	"character",
		  'bam2',	'ba2',	1,	"character",
		  'fld',	'fld',	1,	"character",
		  'meanflnormal','mfn',	2,	"integer",
		  'sdflnormal',	'sfn',	2,	"integer",	
		  'test.type',	'tt', 	1, 	"character",
		  'numsamp',	'n', 	1, 	"integer",
		  'groups',	'grp', 	2, 	"character",
		  'gene.id',	'g',	2, 	"character",
		  'output',	'o', 	1,	"character",
		  'help',	'h', 	0,      "logical"),
		  	ncol=4,byrow=TRUE)

ret.opts<-getopt(options,args)

if ( !is.null(ret.opts$help) ) {
  cat(getopt(options, usage=TRUE));
  q(status=1);
}

# Assignments
transcript.info <- ret.opts$gtf
output.dir <- ret.opts$output

# bam lists
bam.list1 <- list.files(ret.opts$bam1, pattern = "bam$", full.names=TRUE)
bam.list2 <- list.files(ret.opts$bam2, pattern = "bam$", full.names=TRUE)

# test type
if(is.null(ret.opts$testtype))
{
	test.type <- "SKK"
} else
{
	test.type <- ret.opts$test.type
}

if(length(test.type)>1)
{
	test.type <- unlist(strsplit(ret.opts$test.type, ","))
}


# Main function
if((is.null(ret.opts$fld)) || (ret.opts$fld == "empirical"))
{
        FLD <- "empirical"
		IUTA(bam.list1, bam.list2, transcript.info, rep.info.1 = rep(1, length(bam.list1)), rep.info.2 = rep(1, length(bam.list2)), FLD = FLD, test.type = test.type,
    		output.dir = output.dir, output.na = TRUE, genes.interested = "all")
}else if(ret.opts$fld == "normal")
{
        FLD <- "normal"
		mean.FL.normal <- ret.opts$meanflnormal
        sd.FL.normal <- ret.opts$sdflnormal
        IUTA(bam.list1, bam.list2, transcript.info, rep.info.1 = rep(1, length(bam.list1)), rep.info.2 = rep(1, length(bam.list2)), FLD = FLD, mean.FL.normal = mean.FL.normal, 		     sd.FL.normal = sd.FL.normal, test.type = test.type, output.dir = output.dir, output.na = TRUE, genes.interested = "all")
}


# Estimate output 
estimates <- paste(output.dir,"estimates.txt",sep="/")

# Read the estimates file
data <- read.table(estimates, h = T)
genes <- data[,1]
gene.uni <- unique(genes)

# pie_compare and bar_compare
source('/pie_compare.R')
source('/pie_plot.R')
source('/bar_compare.R')

if(!is.null(ret.opts$gene.id))
{
	numb <- ret.opts$numsamp
	gene.name <- ret.opts$gene.id
	group.name <- ret.opts$groups
	group.name <- unlist(strsplit(ret.opts$groups, ","))

	# pie chart
	pie_compare(gene.name, n1 = numb, geometry = "Euclidean", adjust.weight = 1e-2, output.file = paste("Pieplot_", gene.name, ".pdf", sep = ""), group.name = group.name, estimates)
	
	# bar chart
	bar_compare(gene.name, n1 = numb, output.file = paste("Barplot_", gene.name, ".pdf", sep = ""), group.name = group.name, estimates = estimates)
} else
{
	for (gene.name in gene.uni) {

	numb <- ret.opts$numsamp
    	group.name <- ret.opts$groups
    	group.name <- unlist(strsplit(ret.opts$groups, ","))

	# pie chart
	pie_compare(gene.name, n1 = numb, geometry = "Euclidean", adjust.weight = 1e-2, output.file = paste("Pieplot_", gene.name, ".pdf", sep = ""), group.name = group.name, estimates)
	
	# bar chart
	bar_compare(gene.name, n1 = numb, output.file = paste("Barplot_", gene.name, ".pdf", sep = ""), group.name = group.name, estimates = estimates)

 }
	system("tar -zcvf Pieplots.tar.gz Pieplot_*pdf && rm Pieplot_*pdf") 
	system("tar -zcvf Barplots.tar.gz Barplot_*pdf && rm Barplot_*pdf")
}
