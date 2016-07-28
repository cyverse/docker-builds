data <- read.table("final_summary.txt", sep = "\t", head = T)
data$File.Name <- NULL
  
dataMatrix <- data.matrix(data)
dataMatrix <- ifelse(dataMatrix > 0, 1, 0)
  
plot <- apply(dataMatrix, 2, sum)
plot <- round(((plot/max(plot))*100))

png(filename="barplot_lincRNA.png",width=min(3600,1800+800*ncol(plot)/10),height=1800,res=300)
bargraph <- barplot(plot, 
        			main = "LincRNA's across species",
        			col = "Blue",
        			ylab = "Percent A.thaliana homologous lincRNA loci identified",
        			ylim = c(0, max(plot)*1.15),
        			las = 2)
text(bargraph,plot,labels = plot,pos=3,cex=.8)
dev.off()
