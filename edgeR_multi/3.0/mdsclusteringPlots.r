#' Clustering plot
#'
#' Description plots of the counts according to the biological condition
#'
#' @param counts \code{matrix} of counts
#' @param group factor vector of the condition from which each sample belongs
#' @param col colors for the plots (one per biological condition)
#' @return PNG files in the "figures" directory and the matrix of the most expressed sequences
#' @author Upendra Devisetty

mdsclusteringPlots <- function(group=target[,varInt], gene.selection, output.file=output.file, col) {
	
	# Cluster plot
	source("../clusterPlot.R")
	clusterPlot(group=target[,varInt], output.file="cluster.png")  

	# MDS plot
	source("../MDSPlot.R")
	MDSPlot(group=target[,varInt], gene.selection, col, output.file="MDS.png")

	# Heatmap
	source("../heatmap.R")
	Heatmap(output.file="heatmap.png")

	# return(clustplots)

}