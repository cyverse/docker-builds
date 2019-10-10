Heatmap <- function(object=out.DESeq2$dds, output.file="heatmap.png") 
{
                object=out.DESeq2$dds
                counts <- counts(object,normalized=TRUE)
 		counts <- removeNull(counts)

		select <- order(rowMeans(counts),decreasing=TRUE)[1:20]
                
 		nt <- normTransform(object) # defaults to log2(x+1)
		log2.norm.counts <- assay(nt)[select,]
		png(filename=output.file,width=min(3600,1800+800*ncol(counts)/10),height=1500,res=300)
		pheatmap(log2.norm.counts, main="Heatmap")
		dev.off()
}
