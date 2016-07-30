Heatmap <- function(output.file="heatmap.png") 
{
		object=out.edgeR$dge

		# normalization
		dge <- calcNormFactors(object, method=normalizationMethod)
		
		# countspermillion 
		countspermi <- cpm(dge, normalized.lib.sizes=TRUE)

		# Now pick the genes with the top variance over all samples:
		rv <- rowVars(countspermi)
		idx <- order(-rv)[1:20]

		# Plotting
  		png(filename=output.file,width=1800,height=1800,res=300)
		heatmap(countspermi[idx,], Colv = NA)
		title("Heatmap", line = 3)
		dev.off()
}
