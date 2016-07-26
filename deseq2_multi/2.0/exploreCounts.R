#' Explore counts structure
#'
#' Explore counts structure: PCA (DESeq2) or MDS (edgeR) and clustering
#'
#' @param object a \code{DESeqDataSet} from DESeq2 or \code{DGEList} object from edgeR
#' @param group factor vector of the condition from which each sample belongs
#' @param typeTrans transformation method for PCA/clustering with DESeq2: \code{"VST"} or \code{"rlog"}
#' @param gene.selection selection of the features in MDSPlot (\code{"pairwise"} by default)
#' @param col colors used for the PCA/MDS (one per biological condition)
#' @return A list containing the dds object and the results object
#' @author Hugo Varet

exploreCounts <- function(object, group=target[,varInt], typeTrans, gene.selection, col){
    if (typeTrans == "VST") counts.trans <- assay(varianceStabilizingTransformation(object))
    else counts.trans <- assay(rlogTransformation(object))
    source("/PCAPlot.R")
    source("/clusterPlot.R")
    PCAPlot(counts.trans=counts.trans, group=target[,varInt], col=col, output.file="PCAplot.png")
    clusterPlot(counts.trans=counts.trans, group=target[,varInt], output.file="cluster.png")  
  } 
