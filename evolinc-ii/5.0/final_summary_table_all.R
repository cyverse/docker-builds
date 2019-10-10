#!/usr/bin/Rscript
# Script to parse the final summary table to append known lincRNA and/or known gene ID's
# Upendra Kumar Devisetty
# 12/04/16

library(dplyr)


files = list.files(path = "../Homology_Search", pattern = ".csv", full.names = TRUE)
filenames = list.files(path = "../Homology_Search", pattern = "*mod.annotation.*sense.gff", full.names = TRUE)

if (length(files) > 0 && length(filenames) > 0) {
   myfiles = lapply(files, read.delim)  
   data2 = Reduce(function(x, y) merge(x, y, all=TRUE), myfiles)
   full = read.csv("final_summary_table.csv", sep="\t")
   full2 <- full[ , !names(full) %in% names(data2)]
   final_final <- merge(data2, full2, by = 1)
   full.n <- sub("id", "ID", names(full))
   final_final_2 <- final_final[,full.n]
   write.table(final_final_2, file = "final_summary_table.mod.csv", sep = "\t", quote = F, row.names = F)
   
   final <- read.csv("final_summary_table.mod.csv", sep = "\t", stringsAsFactors = F)
   info = file.info(filenames)
   non_empty = rownames(info[info$size != 0, ])  
	 for (f in non_empty) {
    	    data <- read.csv(f, sep = "\t", header = F)
            data2 <- as.data.frame(data %>% group_by(V2) %>% summarise(V18 = sub("^([^[:space:]]+).*","\\1",V18[1])))
    		for (i in 1:nrow(data2)) {
      			sp = unlist(strsplit(as.character(data2[i,1]), "_"))[[1]]
      			query = unlist(strsplit(as.character(data2[i,1]),"_"))[[2]]
      			sub = data2[i,2]
      				if ((sp %in% names(final)) & (query %in% final$ID)) {
        				res1 = final[grep(query, final$ID),][sp][,1]
        				res2 = paste(res1,sub,sep="_")
        				final[grep(query, final$ID),][sp][,1] <- res2
      }
   }
}
  write.table(final, file = "final_summary_table.mod2.csv", quote = F, row.names = F, sep = "\t")	

} else if (length(files) > 0 || length(filenames) < 0) {
  myfiles = lapply(files, read.delim)  
  data2 = Reduce(function(x, y) merge(x, y, all=TRUE), myfiles)
  full = read.csv("final_summary_table.csv", sep="\t")
  full2 <- full[ , !names(full) %in% names(data2)]
  final_final <- merge(data2, full2, by = 1)
  full.n <- sub("id", "ID", names(full))
  final_final_2 <- final_final[,full.n]
  write.table(final_final2, file = "final_summary_table.mod.csv", sep = "\t", quote = F, row.names = F)
  final <- read.csv("final_summary_table.mod.csv", sep = "\t", stringsAsFactors = F)
  write.table(final, file = "final_summary_table.mod2.csv", quote = F, row.names = F, sep = "\t")	

} else if (length(files) < 0 || length(filenames) > 0) {
  final <- read.csv("final_summary_table.csv", sep = "\t", stringsAsFactors = F)
  info = file.info(filenames)
  non_empty = rownames(info[info$size != 0, ])  
  for (f in non_empty) {
    data <- read.csv(f, sep = "\t", header = F)
    data2 <- as.data.frame(data %>% group_by(V2) %>% summarise(V18 = sub("^([^[:space:]]+).*","\\1",V18[1])))
    for (i in 1:nrow(data2)) {
      sp = unlist(strsplit(as.character(data2[i,1]), "_"))[[1]]
      query = unlist(strsplit(as.character(data2[i,1]),"_"))[[2]]
      sub = data2[i,2]
      if ((sp %in% names(final)) & (query %in% final$ID)) {
        res1 = final[grep(query, final$id),][sp][,1]
        res2 = paste(res1,sub,sep="_")
        final[grep(query, final$id),][sp][,1] <- res2
      }
    }
  }
  write.table(final, file = "final_summary_table.mod.csv", quote = F, row.names = F, sep = "\t")
}

