# library
library(dplyr)
library(taxizedb)
library(taxonomizr)
library(seqinr)
library(rentrez)
library(XML)
library(xlsx)
library(ggplot2)
library(Biostrings)
library(rBLAST)

# merge tables
list_of_rds <- list.files(".")
list_of_rds <- list_of_rds[grepl(".rds$", list_of_rds)]
table_list <- list()
table_list_top <- list()
for(i in 1:length(list_of_rds)){
  print(i)
  BLAST_table <- readRDS(list_of_rds[i])
  table_list[[i]] <- BLAST_table
}
merged_table <- do.call(rbind, table_list)
saveRDS(merged_table, file = "Trinity_RNA_sense_proteins_filtered_cas_merged.rds")
