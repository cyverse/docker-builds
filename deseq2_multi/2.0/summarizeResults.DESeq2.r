#' Summarize DESeq2 analysis
#'
#' Summarize DESeq2 analysis: diagnotic plots, dispersions plot, summary of the independent filtering, export results...
#'
#' @param out.DESeq2 the result of \code{run.DESeq2()}
#' @param group factor vector of the condition from which each sample belongs
#' @param independentFiltering \code{TRUE} or \code{FALSE} to perform the independent filtering or not
#' @param cooksCutoff outliers detection threshold (TRUE to let DESeq2 choosing it or FALSE to disable the outliers detection)
#' @param alpha significance threshold to apply to the adjusted p-values
#' @param col colors for the plots
#' @return A list containing: (i) a list of \code{data.frames} from \code{exportResults.DESeq2()}, (ii) the table summarizing the independent filtering procedure and (iii) a table summarizing the number of differentially expressed features
#' @author Hugo Varet

summarizeResults.DESeq2 <- function(out.DESeq2, group=target[,varInt], independentFiltering, cooksCutoff, alpha, col=col){
  
  dds <- out.DESeq2$dds
  results <- out.DESeq2$results
  
  # diagnostic of the size factors
  source("/diagSizeFactorsPlots.R")
  diagSizeFactorsPlots(dds=dds)
  
  # boxplots before and after normalisation
  source("/countsBoxplots.R")
  countsBoxplots(dds, group=group, col=col)
  
  # dispersions plot
  source("/dispersionsPlot.R")
  dispersionsPlot(dds=dds)
  
  # results of the independent filtering
  if (independentFiltering){
    tabIndepFiltering <- tabIndepFiltering(results)
    cat("Number of features discarded by the independent filtering:\n")
    print(tabIndepFiltering, quote=FALSE)
  } else{
    tabIndepFiltering <- NULL
  }
  
  # exporting results of the differential analysis
  source("/exportResults.DESeq2.R")
  complete <- exportResults.DESeq2(out.DESeq2, group=target[,varInt], alpha=alpha)
  
  # small table with number of differentially expressed features
  nDiffTotal <- nDiffTotal(complete=complete, alpha=alpha)
  cat("\nNumber of features down/up and total:\n")
  print(nDiffTotal, quote=FALSE)
  
  # histograms of raw p-values
  source("/rawpHist.R")
  rawpHist(complete=complete, output.file="rawpHist.png")
  
  # MA-plots
  source("/MAPlot.R")
  MAPlot(complete=complete, alpha=alpha, output.file="MAPlot.png")
 
  # Volcano plots
  source("/volcanoPlot.r")
  volcanoPlot(complete=complete, alpha=alpha, output.file="volcanoPlot.png")
 
  return(list(complete=complete, tabIndepFiltering=tabIndepFiltering, nDiffTotal=nDiffTotal))
}
