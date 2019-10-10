#!/usr/bin/Rscript

library(getopt)
library(reshape2)

args<-commandArgs(TRUE)

options<-matrix(c('spe', 's', 1,   "character",
		  'query_sp', 'q', 2,	"character", 	
                  'help', 'h', 0,   "logical"),
                ncol=4,byrow=TRUE)

ret.opts<-getopt(options,args)

if ( !is.null(ret.opts$help) ) {
  cat(getopt(options, usage=TRUE));
  q(status=1);
}

td <- ret.opts$query_sp

td1 <- paste0("^",td,".*$")

species <- ret.opts$spe

species_list <- read.table(species)

myvars <- species_list$V1 %in% td

species_list <- data.frame(species_list$V1[!myvars])

names(species_list) <- "V1"

result = c()

for (f in list.files(path = "lincRNA_families", pattern = ".fasta", full.names = TRUE)) {
  fil <- read.table(f)
  gene.id <- grep(">", fil$V1, value = T)
  gene.id <- sub(">", "", gene.id)
  filter <- grep("TBH_1_", gene.id, value = T)
  filter <- gsub("Known_Gene_Antisense_Known_lincRNA", "KGAKL", filter)
  filter <- gsub("Known_Gene_Sense_Known_lincRNA", "KGSKL", filter)
  filter <- gsub("Known_Gene_Sense$", "KGS", filter)
  filter <- gsub("Known_Gene_Antisense$", "KGA", filter)
  filter <- gsub("Known_lincRNA", "KL", filter)
  result = c(result, filter)
}

res <- c(unique(result))

res <- res[!grepl(td1, res)]

td2 <- paste0(td,"_")

res <- sub(td2, "",res)

new <- matrix(unlist(strsplit(res, "_")), ncol = 5, byrow=T)

sp <- new[,1]
gene <- new[,2]
fun <- new[,5]
  
df1 <- as.data.frame(matrix(NA, ncol = nrow(species_list), nrow = length(gene))) 

for (i in species_list) {colnames(df1) <- i}
t <- data.frame(rep("NA",length(gene)))
colnames(t) <- "id"
df1 <- cbind(t, df1)

df1$id <- gene
df1[cbind(1:nrow(df1), match(sp, names(df1)))] <- fun
df1 <- dcast(subset(melt(df1, id = "id"), value != "Not found"), id ~ variable)
df1[is.na(df1)] <- "Not found"
df1[df1=="KGS"] <- "Known_Gene_Sense"
df1[df1=="KGA"] <- "Known_Gene_Antisense"
df1[df1=="KL"] <- "Known_lincRNA"
df1[df1=="KGSKL"] <- "Known_Gene_Sense_Known_lincRNA"
df1[df1=="KGAKL"] <- "Known_Gene_Antisense_Known_lincRNA"

write.table(df1, "final_summary_table.csv", row.names = F, col.names = T, quote = F, sep = "\t")
