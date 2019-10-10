#!/usr/bin/env Rscript

#.libPaths(getwd())
#install.packages("getopt", repos="http://lib.stat.cmu.edu/R/CRAN", lib=getwd())
#source("http://bioconductor.org/biocLite.R")
#biocLite("ballgown")
#.libPaths(getwd())

library("getopt")
library("ballgown")
library("genefilter")
library("ggplot2")
library("dplyr")
library("reshape2")

#ARGUMENTS
args<-commandArgs(trailingOnly=TRUE)

options<-matrix(c('design', 'd', 1, "character",
                  'covariate', 'o', 1, "character",
                  'datadir', 's',1,"character"),
                   ncol=4, byrow=TRUE)


ret.opts<-getopt(options, args)

designMatrix<-ret.opts$design
pheno_data <- read.table(file=designMatrix, header = TRUE, sep = "\t")
covariateType<-ret.opts$covariate
dataDir<-ret.opts$datadir
sample_full_path<-paste(dataDir,pheno_data$ID, sep = '/')
bg <- ballgown(samples=as.vector(sample_full_path),pData=pheno_data)

# Perform DE with no filtering
results_genes <- stattest(bg, feature='gene', meas='FPKM', covariate=covariateType, getFC=TRUE)
results_genes$logfc <- log2(results_genes$fc)
results_genes <- results_genes[,c(1,2,3,6,4,5)]
results_genes <- arrange(results_genes,pval)

results_transcripts <- stattest(bg, feature='transcript', meas='FPKM', covariate=covariateType, getFC=TRUE)
results_transcripts$logfc <- log2(results_transcripts$fc)
results_transcripts <- results_transcripts[,c(1,2,3,6,4,5)]
results_transcripts <- data.frame(geneNames=geneNames(bg), geneIDs=geneIDs(bg), results_transcripts)
results_transcripts <- arrange(results_transcripts,pval)

write.table(results_transcripts,"results_trans.tsv",sep="\t", row.names=F, quote=F)
write.table(results_genes,"results_gene.tsv",sep="\t", row.names=F, quote=F)

# Plotting
# Histogram plot of transcripts per gene
transcript_gene_table <- indexes(bg)$t2g
counts <- table(transcript_gene_table[,"g_id"])
c_one <- length(which(counts == 1))
c_more_than_one <- length(which(counts > 1))
c_max <- max(counts)
pdf(file = "Histogram_number_transcripts_gene.pdf")
hist(counts, breaks=50, col="bisque4", xlab="Transcripts per gene", main="Distribution of transcript count per gene")
legend_text = c(paste("Genes with one transcript =", c_one), paste("Genes with more than one transcript =", c_more_than_one), paste("Max transcripts for single gene = ", c_max))
legend("topright", legend_text, lty=NULL)
dev.off()

# Histogram of distribution of transcript sizes
full_table <- texpr(bg , 'all')
pdf(file = "Histogram_transcript_sizes.pdf")
hist(full_table$length, breaks=50, xlab="Transcript length (bp)", main="Distribution of transcript lengths", col="steelblue")
dev.off()

# Filter low-abundance genes Here we remove all transcripts with a variance across samples less than one
bg_filt <- subset(bg,"rowVars(texpr(bg)) > 1", genomesubset=TRUE)
results_transcripts <- stattest(bg_filt, feature="transcript", covariate=covariateType, getFC=TRUE, meas="FPKM")
results_transcripts$logfc <- log2(results_transcripts$fc)
results_transcripts <- results_transcripts[,c(1,2,3,6,4,5)]
results_transcripts <- data.frame(geneNames=geneNames(bg_filt), geneIDs=geneIDs(bg_filt), results_transcripts)
results_transcripts <- arrange(results_transcripts,pval)

results_genes <- stattest(bg_filt, feature="gene", covariate=covariateType, getFC=TRUE, meas="FPKM")
results_genes$logfc <- log2(results_genes$fc)
results_genes <- results_genes[,c(1,2,3,6,4,5)]
results_genes <- arrange(results_genes,pval)

write.table(results_transcripts,"results_trans_filter.tsv",sep="\t", row.names=F, quote=F)
write.table(results_genes,"results_gene_filter.tsv",sep="\t", row.names=F, quote=F)

# DE plot
sig <- which(results_genes$qval<0.05)
results_genes[,"de"] <- log2(results_genes[,"fc"])
pdf(file = "DE_plot.pdf")
hist(results_genes[sig,"de"], breaks=50, col="seagreen", xlab="log2(Fold change)", main="Distribution of differential expression values")
abline(v=-2, col="black", lwd=2, lty=2)
abline(v=2, col="black", lwd=2, lty=2)
legend("topleft", "Fold-change > 4", lwd=2, lty=2)
dev.off()

# Identify genes with q value < 0.05
sig_transcripts <- subset(results_transcripts,results_transcripts$qval<0.05)
sig_genes <- subset(results_genes,results_genes$qval<0.05)
write.table(sig_transcripts,"results_trans_filter.sig.tsv",sep="\t", row.names=F, quote=F)
write.table(sig_genes,"results_gene_filter.sig.tsv",sep="\t", row.names=F, quote=F)

#plotting

# MDS plot
fpkm <- gexpr(bg_filt)
short_names <- sub("FPKM.", "", colnames(fpkm))
fpkm <- as.data.frame(fpkm)
data_columns <- 1:ncol(fpkm)
fpkm[,"sum"] <- apply(fpkm[,data_columns], 1, sum)
i <- which(fpkm[,"sum"] > 5)
r <- cor(fpkm[i,data_columns], use="pairwise.complete.obs", method="pearson")
d <- 1-r
mds <- cmdscale(d, k=2, eig=TRUE)
pdf(file = "MDS_plot.pdf")
par(mfrow=c(1,1))
plot(mds$points, type="n", xlab="", ylab="", main="MDS distance plot (all non-zero genes) for all libraries", xlim=c(-0.15,0.15), ylim=c(-0.15,0.15))
points(mds$points[,1], mds$points[,2], col="grey", cex=2, pch=16)
text(mds$points[,1], mds$points[,2], short_names)
dev.off()

# # Boxplot
fpkm <- gexpr(bg_filt)
fpkm_l <- log2(fpkm+1)
fpkm_m <- melt(fpkm_l)
fpkm_m$ID <- gsub("FPKM.", "", fpkm_m$Var2)
head(fpkm_m)
head(pheno_data)
fpkm_m2 <- merge(fpkm_m, pheno_data, by = "ID")
head(fpkm_m2)

p <- ggplot(fpkm_m2, aes(ID, value)) + geom_boxplot() +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(axis.title.x  = element_blank()) +
  labs(y = "log2(FPKM+1)") +
  guides(fill=FALSE) +
  ggtitle("Boxplot of FPKM values")+
  theme(plot.title = element_text(hjust = 0.5))
ggsave(p, file = "FPKM_boxplot.pdf")

# MAplot
results_transcripts$mean <- rowMeans(texpr(bg_filt))

p2 <- ggplot(results_transcripts, aes(log2(mean), log2(fc), colour = qval<0.05)) +
  scale_color_manual(values=c("#999999", "#FF0000")) +
  geom_point() +
  geom_hline(yintercept=0)+
  ggtitle("MA plot") +
  theme(plot.title = element_text(hjust = 0.5))
ggsave(p2, file = "MA_plot.pdf")


