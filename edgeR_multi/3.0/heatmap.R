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
  		png(filename=output.file,width=min(3600,1800+800*ncol(counts)/10),height=1500,res=300)
                pheatmap(countspermi[idx,], main="Heatmap")
		dev.off()
}
