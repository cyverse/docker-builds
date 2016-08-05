#!/usr/bin/Rscript

library(getopt)

args<-commandArgs(TRUE)

options<-matrix(c('spe', 'q', 1,   "character",
                  'help', 'h', 0,   "logical"),
                   ncol=4,byrow=TRUE)

ret.opts<-getopt(options,args)

if ( !is.null(ret.opts$help) ) {
  cat(getopt(options, usage=TRUE));
  q(status=1);
}

species <- ret.opts$spe

dataFile <- read.table("final_summary.txt", sep = "\t", header = TRUE)

dataFile$File.Name <- NULL

dataMatrix <- data.matrix(dataFile)

dataMatrix <- ifelse(dataMatrix > 0, 1, 0)

result <- c()

for (g in colnames(dataMatrix)) { 
  test2 <- strsplit(g, " ")[[1]]
  t1 <- substring(test2, 1, 1)
  t2 <- substring(test2, 2:3)[1]
  t3 <- paste(t1, t2, sep = ".")
  result = append(result,t3)
}

colnames(dataMatrix) <- result

plot <- apply(dataMatrix, 2, sum)

plot <- round(((plot/max(plot))*100))

png(filename="lincRNA_barplot.png",width=min(3600,1800+800*ncol(plot)/10),height=1800,res=300)

string1 = "Percent "
string2 = " homologous lincRNA loci identified"
new = paste0(string1, species, string2) 

bargraph <- barplot(plot, col = "blue",
                    ylab = new, ylim = c(0, max(plot)*1.15),
                    font = 3, yaxt = "n"
                   ,las = 2)

axis(2)
                    
text(bargraph,plot,labels = plot,pos=3,cex=.8)
dev.off()