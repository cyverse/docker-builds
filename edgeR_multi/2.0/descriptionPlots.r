#' Description plots of the counts
#'
#' Description plots of the counts according to the biological condition
#'
#' @param counts \code{matrix} of counts
#' @param group factor vector of the condition from which each sample belongs
#' @param col colors for the plots (one per biological condition)
#' @return PNG files in the "figures" directory and the matrix of the most expressed sequences
#' @author Hugo Varet

descriptionPlots <- function(counts, n=3, group=target[,varInt],output.file=output.file, col){
 # total number of reads per sample
 source("/barplotTotal.R")
 barplotTotal(counts, group=target[,varInt], output.file="barplotTotal.png", col)

 # percentage of null counts per sample
 source("/barplotNull.R")
 barplotNull(counts, group=target[,varInt], output.file="barplotNull.png", col)

 # distribution of counts per sample
 source("/densityPlot.R")
 densityPlot(counts, group=target[,varInt], output.file="densplot.png", col)

 # features which catch the most important number of reads
 source("/majSequences.R")
 majSequences(counts,  n=3, group=target[,varInt], output.file="majSeq.png", col)

 # SERE and pairwise scatter plots
 source("/pairwiseScatterPlots.R")
 cat("Matrix of SERE statistics:\n")
 print(tabSERE(counts))
 pairwiseScatterPlots(counts, group=target[,varInt], output.file="pairwiseScatter.png")
 
 return(majSequences)
}
