#!/usr/bin/Rscript

library(getopt)
library(reshape2)

args<-commandArgs(TRUE)

options<-matrix(c('spe', 's', 1,   "character",
                  'help', 'h', 0,   "logical"),
                   ncol=4,byrow=TRUE)

ret.opts<-getopt(options,args)

if ( !is.null(ret.opts$help) ) {
  cat(getopt(options, usage=TRUE));
  q(status=1);
}

species <- ret.opts$spe

species_list <- read.table(species)

td <- species_list[1,]

species_list <- species_list[-1,]

species_list <- data.frame(species_list)

result = c()

for (f in list.files(path = "lincRNA_families", pattern = ".fasta", full.names = TRUE)) {
  fil <- read.table(f)
  gene.id <- grep(">", fil$V1, value = T)
  gene.id <- sub(">", "", gene.id)
  filter <- grep("[[:print:]]+(_TBH_1)_[[:print:]]", gene.id, value = T)
  result = c(result, filter)
}

res <- c(unique(result))

res <- res[!grepl(td, res)]

sp <- sub("([[:alpha:]]+)_(.*[[:alnum:]]*_TBH_1)_([A-Za-z_]+)","\\1", res)
gene <- sub("([[:alpha:]]+)_(.*[[:alnum:]]*_TBH_1)_([A-Za-z_]+)","\\2", res)
fun <- sub("([[:alpha:]]+)_(.*[[:alnum:]]*_TBH_1)_([A-Za-z_]+)","\\3", res)

df1 <- as.data.frame(matrix(NA, ncol = nrow(species_list), nrow = length(gene))) 
for (i in species_list) {colnames(df1) <- i}
t <- data.frame(rep("NA",length(gene)))
colnames(t) <- "id"
df1 <- cbind(t, df1)

df1$id <- gene
df1[cbind(1:nrow(df1), match(sp, names(df1)))] <- fun
df1 <- dcast(subset(melt(df1, id = "id"), value != "Not found"), id ~ variable + value)
df1[is.na(df1)] <- "Not found"
names(df1) <- sub("([[:alpha:]]+)_([A-Za-z_]+)","\\1", colnames(df1))
write.table(df1, "final_summary_table.tsv", row.names = F, col.names = T, quote = F, sep = "\t")
