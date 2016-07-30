#' Summarize edgeR analysis
#'
#' Summarize edgeR analysis: diagnotic plots, dispersions plot, summary of the independent filtering, export results...
#'
#' @param out.edgeR the result of \code{run.edgeR()}
#' @param group factor vector of the condition from which each sample belongs
#' @param counts matrix of raw counts
#' @param alpha significance threshold to apply to the adjusted p-values
#' @param col colors for the plots
#' @return A list containing: (i) a list of \code{data.frames} from \code{exportResults.edgeR()} and (ii) a table summarizing the number of differentially expressed features
#' @author Hugo Varet

summarizeResults.edgeR <- function(out.edgeR, group=target[,varInt], counts, alpha, col){  
  
  # boxplots before and after normalisation
  source("/countsBoxplots.R")
  countsBoxplots(out.edgeR$dge, col, group=target[,varInt], output.file="countsBoxplots.png")

  # dispersions
  source("/BCVPlot.R")
  BCVPlot(dge=out.edgeR$dge, output.file="BCV.png")
  
  # exporting results of the differential analysis
  source("/exportResults.edgeR.R")
  complete <- exportResults.edgeR(out.edgeR, group=target[,varInt], counts, alpha, OutDir)

  # small table with number of differentially expressed features
  nDiffTotal <- nDiffTotal(complete, alpha=alpha)
  cat("Number of features down/up and total:\n")
  print(nDiffTotal, quote=FALSE)
  
  # histograms of raw p-values
  source("/rawpHist.R")
  rawpHist(complete, output.file="rawpHist.png")
  
  # MA-plots
  source("/MAPlot.R")
  MAPlot(complete, alpha=alpha, output.file="MAPlot.png")
  
  # Volcano plots
  source("/volcanoPlot.r")
  volcanoPlot(complete, alpha=alpha, output.file="volcanoPlot.png")
  
  return(list(complete, nDiffTotal))
}
