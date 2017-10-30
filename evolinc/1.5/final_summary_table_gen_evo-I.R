#!/usr/bin/Rscript

# Install dependencies
library(Biostrings)
library(splitstackshape)
library(dplyr)
library(getopt)
library(splitstackshape)

args<-commandArgs(TRUE)

#############################################################
## Command-line preparations and getopt parsing ##
#############################################################

options<-matrix(c('lincRNA',	  'l',	1,	"character",
		  		        'lincRNAbed',	'b',	1,	"character",
		  		        'overlap',	  'v',	2,	"character",
		  		        'tss',		    't',	2,	"character",
		  		        'help',		    'h', 	0,  "logical"),
		  		         ncol=4, byrow=TRUE)

ret.opts<-getopt(options,args)

if ( !is.null(ret.opts$help) ) {
  cat(getopt(options, usage=TRUE));
  q(status=1);
}


###############################################################################################################
# No TSS or No Overlap
###############################################################################################################

if(is.null(ret.opts$tss) && is.null(ret.opts$overlap)){
  
  # LincRNA fasta file
  fastaFile <- readDNAStringSet(ret.opts$lincRNA)
  lincRNA_ID = names(fastaFile)
  sequence = paste(fastaFile)
  size = width(fastaFile)
  df <- data.frame(lincRNA_ID, size)

  # Bed file
  # split concatenated column by `;`
  bed_File <- read.csv(ret.opts$lincRNAbed, sep = "\t", header = F)

  # Search the tss column
  new <- bed_File[grep("tss_id", bed_File$V10),]

if (nrow(new) == 0) {
  # Without tss
  new2 <- bed_File[!grepl("tss_id", bed_File$V10),]
  results2 <- cSplit(new2, splitCols= "V10", sep = ";", drop = TRUE)
  results2 <- as.data.frame(results2)
  gene2_co <- which(apply(results2, 2, function(x) any(grepl("XLOC_*", x))))
  TCONS2_co <- which(apply(results2, 2, function(x) any(grepl("TCONS*", x))))
  exon2_co <- which(apply(results2, 2, function(x) any(grepl("exon_number", x))))
  cc2_co <- which(apply(results2, 2, function(x) any(grepl("class_code", x))))
  comb2 <- cbind(results2[gene2_co], results2[TCONS2_co], results2[exon2_co], results2[cc2_co])
  comb2.m <- as.data.frame(cSplit(comb2, splitCols = c("V10_1", "V10_5", "V10_3"), sep = " ", drop = TRUE))
  comb2.m <- comb2.m[,c(3,5,7,1)]
  names(comb2.m) <- c("gene", "transcript", "exon_numb", "class_code")
  
  t <- paste0(comb2.m$transcript, ".gene=", comb2.m$gene)
  bind <- as.data.frame(cbind(t,comb2.m$exon_numb))
  bind$t <- as.character(bind$t)
  bind$V2 <- as.character(bind$V2)
  bedfinal <- as.data.frame(bind %>% group_by(t) %>% summarise(V2 = max(V2)))
  
  # non present columns
  cols <- ncol(bedfinal)
  cols <- cols + 1 
  bedfinal[,cols] <- NA
  cols <- cols + 1 
  bedfinal[,cols] <- NA
  cols <- cols + 1 
  bedfinal[,cols] <- NA
  merge1 <- merge(df, bedfinal, by=1)
  merge1 <- merge1[,c(1,2,4:6,3)]
  names(merge1) <- c("lincRNA_ID", "Size(bp)", "Overlapping_known_lincRNA", "gene_ID", "Has_TSS_data", "Number_of_exons") 
  write.table(merge1, file = "final_Summary_table.tsv", row.names = F, col.names = T, quote = F, sep = "\t")
  
} else {
  # With tss
  results <- as.data.frame(cSplit(new, splitCols= "V10", sep = ";", drop = TRUE))
  gene_co <- which(apply(results, 2, function(x) any(grepl("XLOC_*", x))))
  TCONS_co <- which(apply(results, 2, function(x) any(grepl("TCONS*", x))))
  exon_co <- which(apply(results, 2, function(x) any(grepl("exon_number", x))))
  cc_co <- which(apply(results, 2, function(x) any(grepl("class_code", x))))
  comb <- cbind(results[gene_co], results[TCONS_co], results[exon_co], results[cc_co])
  comb.m <- as.data.frame(cSplit(comb, splitCols = c(1, 2, 3), sep = " ", drop = TRUE))
  comb.m <- comb.m[,c(3,5,7,1)]
  names(comb.m) <- c("gene", "transcript", "exon_numb", "class_code")
  
  # Without tss
  new2 <- bed_File[!grepl("tss_id", bed_File$V10),]
  
  if (nrow(new2) != 0) {
    results2 <- cSplit(new2, splitCols= "V10", sep = ";", drop = TRUE)
    results2 <- as.data.frame(results2)
    gene2_co <- which(apply(results2, 2, function(x) any(grepl("XLOC_*", x))))
    TCONS2_co <- which(apply(results2, 2, function(x) any(grepl("TCONS*", x))))
    exon2_co <- which(apply(results2, 2, function(x) any(grepl("exon_number", x))))
    cc2_co <- which(apply(results2, 2, function(x) any(grepl("class_code", x))))
    comb2 <- cbind(results2[gene2_co], results2[TCONS2_co], results2[exon2_co], results2[cc2_co])
    comb2.m <- as.data.frame(cSplit(comb2, splitCols = c(1, 2, 3), sep = " ", drop = TRUE))
    comb2.m <- comb2.m[,c(3,5,7,1)]
    names(comb2.m) <- c("gene", "transcript", "exon_numb", "class_code")
    
    # Combine now
    results_all <- rbind(comb.m, comb2.m)
    head(results_all)
    
    # Adding gene, transcript and exon number columns
    t <- paste0(results_all$transcript, ".gene=", results_all$gene)
    bind <- as.data.frame(cbind(t,results_all$exon_numb))
    bind$t <- as.character(bind$t)
    bind$V2 <- as.character(bind$V2)
    bedfinal <- as.data.frame(bind %>% group_by(t) %>% summarise(V2 = max(V2)))
    
    # non present columns
    cols <- ncol(bedfinal)
    cols <- cols + 1 
    bedfinal[,cols] <- NA
    cols <- cols + 1 
    bedfinal[,cols] <- NA
    cols <- cols + 1 
    bedfinal[,cols] <- NA
    merge1 <- merge(df, bedfinal, by=1)
    merge1 <- merge1[,c(1,2,4:6,3)]
    names(merge1) <- c("lincRNA_ID", "Size(bp)", "Overlapping_known_lincRNA", "gene_ID", "Has_TSS_data", "Number_of_exons") 
    write.table(merge1, file = "final_Summary_table.tsv", row.names = F, col.names = T, quote = F, sep = "\t")
 
  } else {
    t <- paste0(comb.m$transcript, ".gene=", comb.m$gene)
    bind <- as.data.frame(cbind(t,comb.m$exon_numb))
    bind$t <- as.character(bind$t)
    bind$V2 <- as.character(bind$V2)
    bedfinal <- as.data.frame(bind %>% group_by(t) %>% summarise(V2 = max(V2)))
    
    # non present columns
    cols <- ncol(bedfinal)
    cols <- cols + 1 
    bedfinal[,cols] <- NA
    cols <- cols + 1 
    bedfinal[,cols] <- NA
    cols <- cols + 1 
    bedfinal[,cols] <- NA
    merge1 <- merge(df, bedfinal, by=1)
    merge1 <- merge1[,c(1,2,4:6,3)]
    names(merge1) <- c("lincRNA_ID", "Size(bp)", "Overlapping_known_lincRNA", "gene_ID", "Has_TSS_data", "Number_of_exons")
    write.table(merge1, file = "final_Summary_table.tsv", row.names = F, col.names = T, quote = F, sep = "\t")
  }       
}

###############################################################################################################
# No TSS file
###############################################################################################################

} else if(is.null(ret.opts$tss)){
  # Fasta file
  fastaFile <- readDNAStringSet(ret.opts$lincRNA)
  lincRNA_ID = names(fastaFile)
  sequence = paste(fastaFile)
  size = width(fastaFile)
  df <- data.frame(lincRNA_ID, size)
  
  # Overlapping lincRNA's
  fastaFile_over <- readDNAStringSet(ret.opts$overlap)
  lincRNA_ID_over <- names(fastaFile_over)
  len <- length(lincRNA_ID_over)
  newmat <- matrix(lincRNA_ID_over, ncol =1 , nrow = len, byrow = T)
  overlap <- cSplit(as.data.table(newmat), "V1", "_o")

    # Finding out the lincRNA gene id:
    info = file.info("./intersect_output2.txt")
    if (info$size == 0) {
      data4 <- data.frame(gene_id = rep("NA",nrow(df)))
      }else {    
        data1 <- read.table("./intersect_output2.txt", sep="\t")
        names(data1) <- c("ID","gene")
        data1$gene <- as.character(data1$gene)
        gene2 <- unlist(strsplit(data1$gene, "_"))
        data2 <- matrix(gene2, ncol = 4, byrow = T)
        data3 <- as.data.frame(cbind(data1[1], data2[,2]))
        names(data3) <- c("id", "gene")
        data4 <- as.data.frame(data3 %>% group_by(id) %>% summarise(gene[1]))
        names(data4)[2] <- "gene_id" 
    }

    # Bed file
    # split concatenated column by `;`
    bed_File <- read.csv("lincRNAs.bed", sep = "\t", header = F)

    # Search the tss column
    new <- bed_File[grep("tss_id", bed_File$V10),]

    if (nrow(new) == 0) {
      # Without tss
      new2 <- bed_File[!grepl("tss_id", bed_File$V10),]
      results2 <- cSplit(new2, splitCols= "V10", sep = ";", drop = TRUE)
      results2 <- as.data.frame(results2)
      gene2_co <- which(apply(results2, 2, function(x) any(grepl("XLOC_*", x))))
      TCONS2_co <- which(apply(results2, 2, function(x) any(grepl("TCONS*", x))))
      exon2_co <- which(apply(results2, 2, function(x) any(grepl("exon_number", x))))
      cc2_co <- which(apply(results2, 2, function(x) any(grepl("class_code", x))))
      comb2 <- cbind(results2[gene2_co], results2[TCONS2_co], results2[exon2_co], results2[cc2_co])
      comb2.m <- as.data.frame(cSplit(comb2, splitCols = c("V10_1", "V10_5", "V10_3"), sep = " ", drop = TRUE))
      comb2.m <- comb2.m[,c(3,5,7,1)]
      names(comb2.m) <- c("gene", "transcript", "exon_numb", "class_code")
      
      t <- paste0(comb2.m$transcript, ".gene=", comb2.m$gene)
      bind <- as.data.frame(cbind(t,comb2.m$exon_numb))
      bind$t <- as.character(bind$t)
      bind$V2 <- as.character(bind$V2)
      bedfinal <- as.data.frame(bind %>% group_by(t) %>% summarise(V2 = max(V2)))
      
      # non present columns
      cols <- ncol(bedfinal)
      cols <- cols + 1 
      bedfinal[,cols] <- NA

    } else {
      # With tss
      results <- as.data.frame(cSplit(new, splitCols= "V10", sep = ";", drop = TRUE))
      gene_co <- which(apply(results, 2, function(x) any(grepl("XLOC_*", x))))
      TCONS_co <- which(apply(results, 2, function(x) any(grepl("TCONS*", x))))
      exon_co <- which(apply(results, 2, function(x) any(grepl("exon_number", x))))
      cc_co <- which(apply(results, 2, function(x) any(grepl("class_code", x))))
      comb <- cbind(results[gene_co], results[TCONS_co], results[exon_co], results[cc_co])
      comb.m <- as.data.frame(cSplit(comb, splitCols = c(1, 2, 3), sep = " ", drop = TRUE))
      comb.m <- comb.m[,c(3,5,7,1)]
      names(comb.m) <- c("gene", "transcript", "exon_numb", "class_code")
      
      # Without tss
      new2 <- bed_File[!grepl("tss_id", bed_File$V10),]
      
      if (nrow(new2) != 0) {
        results2 <- cSplit(new2, splitCols= "V10", sep = ";", drop = TRUE)
        results2 <- as.data.frame(results2)
        gene2_co <- which(apply(results2, 2, function(x) any(grepl("XLOC_*", x))))
        TCONS2_co <- which(apply(results2, 2, function(x) any(grepl("TCONS*", x))))
        exon2_co <- which(apply(results2, 2, function(x) any(grepl("exon_number", x))))
        cc2_co <- which(apply(results2, 2, function(x) any(grepl("class_code", x))))
        comb2 <- cbind(results2[gene2_co], results2[TCONS2_co], results2[exon2_co], results2[cc2_co])
        comb2.m <- as.data.frame(cSplit(comb2, splitCols = c(1, 2, 3), sep = " ", drop = TRUE))
        comb2.m <- comb2.m[,c(3,5,7,1)]
        names(comb2.m) <- c("gene", "transcript", "exon_numb", "class_code")
        
        # Combine now
        results_all <- rbind(comb.m, comb2.m)
        head(results_all)
        
        # Adding gene, transcript and exon number columns
        t <- paste0(results_all$transcript, ".gene=", results_all$gene)
        bind <- as.data.frame(cbind(t,results_all$exon_numb))
        bind$t <- as.character(bind$t)
        bind$V2 <- as.character(bind$V2)
        bedfinal <- as.data.frame(bind %>% group_by(t) %>% summarise(V2 = max(V2)))
        
        # non present columns
        cols <- ncol(bedfinal)
        cols <- cols + 1 
        bedfinal[,cols] <- NA

      } else {
        t <- paste0(comb.m$transcript, ".gene=", comb.m$gene)
        bind <- as.data.frame(cbind(t,comb.m$exon_numb))
        bind$t <- as.character(bind$t)
        bind$V2 <- as.character(bind$V2)
        bedfinal <- as.data.frame(bind %>% group_by(t) %>% summarise(V2 = max(V2)))
        
        # non present columns
        cols <- ncol(bedfinal)
        cols <- cols + 1 
        bedfinal[,cols] <- NA
      }       
    }

    merge1 <- merge(x = df, y = overlap, by = 1, all = TRUE)
    if (info$size == 0) {
      merge2 <- cbind(merge1, data4) 
      } else {
      merge2 <- merge(x = merge1, y = data4, by = 1, all = TRUE)
    }

    merge2$V1_2 <- as.character(merge2$V1_2)
    merge2$gene_id <- as.character(merge2$gene_id)
    merge2$V1_2[!is.na(merge2$V1_2)] <- "Yes"
    merge2$V1_2[is.na(merge2$V1_2)] <- "No" 
    merge2$gene_id[is.na(merge2$gene_id)] <- "NA"       
    colnames(merge2)[2] <- "Size(bp)"
    colnames(merge2)[3] <- "Overlapping_known_lincRNA"
    merge3 <- merge(x = merge2, y = bedfinal, by = 1, all = TRUE)
    colnames(merge3)[5] <- "Number_of_exons"
    colnames(merge3)[6] <- "Has_TSS_data"
    merge3 <- merge3[,c(1:4,6,5)]
    # Writing data
    write.table(merge3, file = "final_Summary_table.tsv", row.names = F, col.names = T, quote = F, sep = "\t")


###############################################################################################################
# No Overlap file
###############################################################################################################

} else if(is.null(ret.opts$overlap)) {
    # Fasta file
    fastaFile <- readDNAStringSet(ret.opts$lincRNA)
    lincRNA_ID = names(fastaFile)
    sequence = paste(fastaFile)
    size = width(fastaFile)
    df <- data.frame(lincRNA_ID, size)

    # CAGE supported lincRNA's
    fastaFile_Cage <- readDNAStringSet(ret.opts$tss)
    lincRNA_ID_Cage <- names(fastaFile_Cage)
    len1 <- length(lincRNA_ID_Cage)
    newmat1 <- matrix(lincRNA_ID_Cage, ncol =1 , nrow = len1, byrow = T)
    cage <- cSplit(as.data.table(newmat1), "V1", "_C")

    # Bed file
    # split concatenated column by `;`
    bed_File <- read.csv(ret.opts$lincRNAbed, sep = "\t", header = F)

    # Search the tss column
    new <- bed_File[grep("tss_id", bed_File$V10),]

    if (nrow(new) == 0) {
      # Without tss
      new2 <- bed_File[!grepl("tss_id", bed_File$V10),]
      results2 <- cSplit(new2, splitCols= "V10", sep = ";", drop = TRUE)
      results2 <- as.data.frame(results2)
      gene2_co <- which(apply(results2, 2, function(x) any(grepl("XLOC_*", x))))
      TCONS2_co <- which(apply(results2, 2, function(x) any(grepl("TCONS*", x))))
      exon2_co <- which(apply(results2, 2, function(x) any(grepl("exon_number", x))))
      cc2_co <- which(apply(results2, 2, function(x) any(grepl("class_code", x))))
      comb2 <- cbind(results2[gene2_co], results2[TCONS2_co], results2[exon2_co], results2[cc2_co])
      comb2.m <- as.data.frame(cSplit(comb2, splitCols = c("V10_1", "V10_5", "V10_3"), sep = " ", drop = TRUE))
      comb2.m <- comb2.m[,c(3,5,7,1)]
      names(comb2.m) <- c("gene", "transcript", "exon_numb", "class_code")
      
      t <- paste0(comb2.m$transcript, ".gene=", comb2.m$gene)
      bind <- as.data.frame(cbind(t,comb2.m$exon_numb))
      bind$t <- as.character(bind$t)
      bind$V2 <- as.character(bind$V2)
      bedfinal <- as.data.frame(bind %>% group_by(t) %>% summarise(V2 = max(V2)))
      
      # non present columns
      cols <- ncol(bedfinal)
      cols <- cols + 1 
      bedfinal[,cols] <- NA
      cols <- cols + 1 
      bedfinal[,cols] <- NA

    } else {
      # With tss
      results <- as.data.frame(cSplit(new, splitCols= "V10", sep = ";", drop = TRUE))
      gene_co <- which(apply(results, 2, function(x) any(grepl("XLOC_*", x))))
      TCONS_co <- which(apply(results, 2, function(x) any(grepl("TCONS*", x))))
      exon_co <- which(apply(results, 2, function(x) any(grepl("exon_number", x))))
      cc_co <- which(apply(results, 2, function(x) any(grepl("class_code", x))))
      comb <- cbind(results[gene_co], results[TCONS_co], results[exon_co], results[cc_co])
      comb.m <- as.data.frame(cSplit(comb, splitCols = c(1, 2, 3), sep = " ", drop = TRUE))
      comb.m <- comb.m[,c(3,5,7,1)]
      names(comb.m) <- c("gene", "transcript", "exon_numb", "class_code")
      
      # Without tss
      new2 <- bed_File[!grepl("tss_id", bed_File$V10),]
      
      if (nrow(new2) != 0) {
        results2 <- cSplit(new2, splitCols= "V10", sep = ";", drop = TRUE)
        results2 <- as.data.frame(results2)
        gene2_co <- which(apply(results2, 2, function(x) any(grepl("XLOC_*", x))))
        TCONS2_co <- which(apply(results2, 2, function(x) any(grepl("TCONS*", x))))
        exon2_co <- which(apply(results2, 2, function(x) any(grepl("exon_number", x))))
        cc2_co <- which(apply(results2, 2, function(x) any(grepl("class_code", x))))
        comb2 <- cbind(results2[gene2_co], results2[TCONS2_co], results2[exon2_co], results2[cc2_co])
        comb2.m <- as.data.frame(cSplit(comb2, splitCols = c(1, 2, 3), sep = " ", drop = TRUE))
        comb2.m <- comb2.m[,c(3,5,7,1)]
        names(comb2.m) <- c("gene", "transcript", "exon_numb", "class_code")
        
        # Combine now
        results_all <- rbind(comb.m, comb2.m)
        head(results_all)
        
        # Adding gene, transcript and exon number columns
        t <- paste0(results_all$transcript, ".gene=", results_all$gene)
        bind <- as.data.frame(cbind(t,results_all$exon_numb))
        bind$t <- as.character(bind$t)
        bind$V2 <- as.character(bind$V2)
        bedfinal <- as.data.frame(bind %>% group_by(t) %>% summarise(V2 = max(V2)))
        
        # non present columns
        cols <- ncol(bedfinal)
        cols <- cols + 1 
        bedfinal[,cols] <- NA
        cols <- cols + 1 
        bedfinal[,cols] <- NA

      } else {
        t <- paste0(comb.m$transcript, ".gene=", comb.m$gene)
        bind <- as.data.frame(cbind(t,comb.m$exon_numb))
        bind$t <- as.character(bind$t)
        bind$V2 <- as.character(bind$V2)
        bedfinal <- as.data.frame(bind %>% group_by(t) %>% summarise(V2 = max(V2)))
        
        # non present columns
        cols <- ncol(bedfinal)
        cols <- cols + 1 
        bedfinal[,cols] <- NA
        cols <- cols + 1 
        bedfinal[,cols] <- NA
      }       
    }

    # Merging All three files
    merge1 <- merge(x = df, y = cage, by = 1, all = TRUE)
    merge2 <- merge(x = merge1, y = bedfinal, by = 1, all = TRUE)
    merge2$V1_2 <- as.character(merge2$V1_2)
    merge2$V1_2[is.na(merge2$V1_2)] <- "No"
    merge2$V1_2 <- sub("AGE_PLUS", "Yes", merge2$V1_2)
    merge2$V1_2 <- as.factor(merge2$V1_2)
    colnames(merge2)[2] <- "Size(bp)"
    colnames(merge2)[3] <- "Has_TSS_data"
    colnames(merge2)[4] <- "Number_of_exons"
    colnames(merge2)[5] <- "Overlapping_known_lincRNA"
    colnames(merge2)[6] <- "gene_id"
    # merge2[,order(colnames(merge2))]
    merge2 <- merge2[,c(1:2,5:6,3:4)]
    # Writing data
    write.table(merge2, file = "final_Summary_table.tsv", row.names = F, col.names = T, quote = F, sep = "\t")

########################################################################################################################  
# All present
########################################################################################################################

} else {
    # All lincRNA's
    fastaFile <- readDNAStringSet(ret.opts$lincRNA)
    lincRNA_ID = names(fastaFile)
    sequence = paste(fastaFile)
    size = width(fastaFile)
    df <- data.frame(lincRNA_ID, size)
    
    # Overlapping lincRNA's
    fastaFile_over <- readDNAStringSet(ret.opts$overlap)
    lincRNA_ID_over <- names(fastaFile_over)
    len <- length(lincRNA_ID_over)
    newmat <- matrix(lincRNA_ID_over, ncol =1 , nrow = len, byrow = T)
    overlap <- cSplit(as.data.table(newmat), "V1", "_o")
    
    # CAGE supported lincRNA's
    fastaFile_Cage <- readDNAStringSet(ret.opts$tss)
    lincRNA_ID_Cage <- names(fastaFile_Cage)
    len1 <- length(lincRNA_ID_Cage)
    newmat1 <- matrix(lincRNA_ID_Cage, ncol =1 , nrow = len1, byrow = T)
    cage <- cSplit(as.data.table(newmat1), "V1", "_C")
    
    # Finding out the lincRNA gene id:
    info = file.info("./intersect_output2.txt")
    if (info$size == 0) {
      data4 <- data.frame(gene_id = rep("NA",nrow(df)))
      } else {    
          data1 <- read.table("./intersect_output2.txt", sep="\t")
          names(data1) <- c("ID","gene")
          data1$gene <- as.character(data1$gene)
          gene2 <- unlist(strsplit(data1$gene, "_"))
          data2 <- matrix(gene2, ncol = 4, byrow = T)
          data3 <- as.data.frame(cbind(data1[1], data2[,2]))
          names(data3) <- c("id", "gene")
          data4 <- as.data.frame(data3 %>% group_by(id) %>% summarise(gene[1]))
          names(data4)[2] <- "gene_id" 
    }   
    # Bed file
    # split concatenated column by `;`
    bed_File <- read.csv(ret.opts$lincRNAbed, sep = "\t", header = F)

    # Search the tss column
    new <- bed_File[grep("tss_id", bed_File$V10),]

    if (nrow(new) == 0) {
      # Without tss
      new2 <- bed_File[!grepl("tss_id", bed_File$V10),]
      results2 <- cSplit(new2, splitCols= "V10", sep = ";", drop = TRUE)
      results2 <- as.data.frame(results2)
      gene2_co <- which(apply(results2, 2, function(x) any(grepl("XLOC_*", x))))
      TCONS2_co <- which(apply(results2, 2, function(x) any(grepl("TCONS*", x))))
      exon2_co <- which(apply(results2, 2, function(x) any(grepl("exon_number", x))))
      cc2_co <- which(apply(results2, 2, function(x) any(grepl("class_code", x))))
      comb2 <- cbind(results2[gene2_co], results2[TCONS2_co], results2[exon2_co], results2[cc2_co])
      comb2.m <- as.data.frame(cSplit(comb2, splitCols = c("V10_1", "V10_5", "V10_3"), sep = " ", drop = TRUE))
      comb2.m <- comb2.m[,c(3,5,7,1)]
      names(comb2.m) <- c("gene", "transcript", "exon_numb", "class_code")
      
      t <- paste0(comb2.m$transcript, ".gene=", comb2.m$gene)
      bind <- as.data.frame(cbind(t,comb2.m$exon_numb))
      bind$t <- as.character(bind$t)
      bind$V2 <- as.character(bind$V2)
      bedfinal <- as.data.frame(bind %>% group_by(t) %>% summarise(V2 = max(V2)))

    } else {
      # With tss
      results <- as.data.frame(cSplit(new, splitCols= "V10", sep = ";", drop = TRUE))
      gene_co <- which(apply(results, 2, function(x) any(grepl("XLOC_*", x))))
      TCONS_co <- which(apply(results, 2, function(x) any(grepl("TCONS*", x))))
      exon_co <- which(apply(results, 2, function(x) any(grepl("exon_number", x))))
      cc_co <- which(apply(results, 2, function(x) any(grepl("class_code", x))))
      comb <- cbind(results[gene_co], results[TCONS_co], results[exon_co], results[cc_co])
      comb.m <- as.data.frame(cSplit(comb, splitCols = c(1, 2, 3), sep = " ", drop = TRUE))
      comb.m <- comb.m[,c(3,5,7,1)]
      names(comb.m) <- c("gene", "transcript", "exon_numb", "class_code")
      
      # Without tss
      new2 <- bed_File[!grepl("tss_id", bed_File$V10),]
      
      if (nrow(new2) != 0) {
        results2 <- cSplit(new2, splitCols= "V10", sep = ";", drop = TRUE)
        results2 <- as.data.frame(results2)
        gene2_co <- which(apply(results2, 2, function(x) any(grepl("XLOC_*", x))))
        TCONS2_co <- which(apply(results2, 2, function(x) any(grepl("TCONS*", x))))
        exon2_co <- which(apply(results2, 2, function(x) any(grepl("exon_number", x))))
        cc2_co <- which(apply(results2, 2, function(x) any(grepl("class_code", x))))
        comb2 <- cbind(results2[gene2_co], results2[TCONS2_co], results2[exon2_co], results2[cc2_co])
        comb2.m <- as.data.frame(cSplit(comb2, splitCols = c(1, 2, 3), sep = " ", drop = TRUE))
        comb2.m <- comb2.m[,c(3,5,7,1)]
        names(comb2.m) <- c("gene", "transcript", "exon_numb", "class_code")
        
        # Combine now
        results_all <- rbind(comb.m, comb2.m)
        head(results_all)
        
        # Adding gene, transcript and exon number columns
        t <- paste0(results_all$transcript, ".gene=", results_all$gene)
        bind <- as.data.frame(cbind(t,results_all$exon_numb))
        bind$t <- as.character(bind$t)
        bind$V2 <- as.character(bind$V2)
        bedfinal <- as.data.frame(bind %>% group_by(t) %>% summarise(V2 = max(V2)))
        
      } else {
        t <- paste0(comb.m$transcript, ".gene=", comb.m$gene)
        bind <- as.data.frame(cbind(t,comb.m$exon_numb))
        bind$t <- as.character(bind$t)
        bind$V2 <- as.character(bind$V2)
        bedfinal <- as.data.frame(bind %>% group_by(t) %>% summarise(V2 = max(V2)))
        
      }       
    }
    
    # Merging All three files
    merge1 <- merge(x = df, y = overlap, by = 1, all = TRUE)
    if (info$size == 0) {
        merge2 <- cbind(merge1, data4)}
    else {
        merge2 <- merge(x = merge1, y = data4, by = 1, all = TRUE)
    }

    merge3 <- merge(merge2, cage, by=1, all = TRUE)
    merge4 <- merge(merge3, bedfinal, by=1, all = TRUE)
    names(merge4) <- c("lincRNA_ID", "Size(bp)", "Overlapping_known_lincRNA", "gene_ID", "Has_TSS_data", "Number_of_exons") 
    merge4$gene_ID <- as.character(merge4$gene_ID)
    merge4$gene_ID[is.na(merge4$gene_ID)] <- "NA"
    merge4$Overlapping_known_lincRNA <- as.character(merge4$Overlapping_known_lincRNA)
    merge4$Has_TSS_data <- as.character(merge4$Has_TSS_data)
    merge4$Overlapping_known_lincRNA[!is.na(merge4$Overlapping_known_lincRNA)] <- "Yes"
    merge4$Has_TSS_data[!is.na(merge4$Has_TSS_data)] <- "Yes"
    merge4$Overlapping_known_lincRNA[is.na(merge4$Overlapping_known_lincRNA)] <- "No"
    merge4$Has_TSS_data[is.na(merge4$Has_TSS_data)] <- "No"

    # Writing data
    write.table(merge4, file = "final_Summary_table.tsv", row.names = F, col.names = T, quote = F, sep = "\t")
}