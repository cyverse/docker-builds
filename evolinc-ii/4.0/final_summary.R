#!/usr/bin/Rscript

library(getopt)

args<-commandArgs(TRUE)

options<-matrix(c('spe', 's', 1,   "character",
                  'query_sp', 'q', 1, "character",
                  'help', 'h', 0,   "logical"),
                   ncol=4,byrow=TRUE)

ret.opts<-getopt(options,args)

if ( !is.null(ret.opts$help) ) {
  cat(getopt(options, usage=TRUE));
  q(status=1);
}

species <- read.table(ret.opts$spe)

a1 <- substring(ret.opts$query_sp,1,1)

a2 <- substring(ret.opts$query_sp,2:3)[1]

a3 <- paste(a1, a2, sep = ".")

dataFile <- read.table("final_summary.txt", sep = "\t", header = TRUE)

dataFile$File.Name <- NULL

query_column <- dataFile[ret.opts$query_sp]

dataFile[ret.opts$query_sp] <- NULL

dataFile <- cbind(query_column, dataFile)

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

res1 = paste0("(", plot, ")")

plot1 <- round(((plot/max(plot))*100))

res2 = paste(plot1, res1, sep = " ")

png(filename="lincRNA_barplot.png",width=min(3600,1800+800*ncol(plot)/10),height=1800,res=300)

string1 = "Percent "
string2 = " homologous lincRNA loci identified"
new = paste0(string1, a3, string2) 

bargraph <- barplot(plot1, col = "blue",
                    ylab = new, ylim = c(0, max(plot1)*1.15),
                    font = 3, yaxt = "n"
                   ,las = 2)

axis(2)
                    
text(bargraph, plot1, labels=res2, pos=3, cex=.8)
dev.off()
