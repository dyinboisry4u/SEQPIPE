library(readr)
library(knitr)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(readr)
library(tximport)
library(tibble)
args=commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
  stop("You must provide three arguments: dataDir, quantMethod , ENSG2T_anno (gencode.v39.ENSG2T.annotation.csv).")
}
dataDir=args[1]
quantMethod=args[2]
ENSG2T_anno=args[3]
##merge *.count *.sf 
if (quantMethod == 'featureCounts') {
  samples = read.table(file.path(dataDir,"sampleinfo.txt"),header = F)
  files = file.path(dataDir, "04_quantification", paste0(samples$V2,"_featureCounts.counts"))
  names(files) <- samples$V2
  
  data_list <- list()
  for (file_label in names(files)) {
    file_path <- files[[file_label]]
    
    df <- read.table(file_path, header = T, sep = "\t", colClasses = c("character", "NULL", "NULL", "NULL", "NULL", "NULL", "numeric"))
    colnames(df) <- c("Gene", file_label)
    
    data_list[[file_label]] <- df
  }
  merged_data <- Reduce(function(x, y) merge(x, y, by = "Gene", all = TRUE), data_list)
  write.table(merged_data, file = paste0(dataDir, "/04_quantification/merge.count"), sep = "\t", quote = FALSE, row.names = F,col.names = T)
  
} else if (quantMethod == 'Salmon') {
  
  tx2gene <- read_csv(file.path(ENSG2T_anno))
  samples <- read.table(file.path(dataDir, "sampleinfo.txt"), header = FALSE)
  files <- file.path(dataDir, "04_quantification", paste0(samples$V2, "_salmon_quant"), "quant.sf")
  names(files) <- samples$V2
  txi <- tximport(files, type = "salmon", tx2gene = tx2gene)
  final_counts <- as.data.frame(txi$counts) 
  final_counts=rownames_to_column(final_counts)
  write.table(final_counts, file = paste0(dataDir, "/04_quantification/txim_merge.count"), sep = "\t", quote = FALSE, row.names = F)
} else {
  message("Invalid quantMethod type!")
}

