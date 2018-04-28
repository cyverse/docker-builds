#!/usr/bin/Rscript

library(getopt)

args <- commandArgs(TRUE)

options <- matrix(c('Dir', 'r', 2, "character",
					'Outfile', 'w', 1, "character",
					'help',   'h',    0,      "logical"),
                     ncol=4,byrow=TRUE)

# Usage
# Rscript merge_htseq.R --Dir test --Outfile output_test.txt

ret.opts<-getopt(options,args)

if ( !is.null(ret.opts$help) ) {
  cat(getopt(options, usage=TRUE));
  q(status=1);
}

rawDir <- ret.opts$Dir
Outfile <- ret.opts$Outfile

for (i in list.files(path = rawDir)) {
  final2 <- do.call(rbind, lapply(list.files(path = rawDir, pattern = i, full.names = T), read.table))
  names(final2) <- c("gene", i)
  new_final <- paste("comb_test", i, sep = "_")
  write.table(final2, file = new_final, quote = F, row.names = F, col.names = T, sep = "\t")
  }

df_list <- list()

for (i in list.files(path = ".", pattern = "comb_test", full.names = T)) {
  df_list[[i]] <- read.csv(i, sep = "\t")
}

merged_df <- Reduce(function(x, y) merge(x, y, by = "gene", all = TRUE), df_list)

write.table(merged_df, file = Outfile, quote = F, row.names = F, sep = "\t")

system("rm comb_test*.txt")