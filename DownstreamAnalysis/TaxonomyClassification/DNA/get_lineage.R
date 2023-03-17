# library
library(dplyr)
library(ComplexHeatmap)
library(ggplot2)
library(reshape)
library(taxizedb)

###
# get lineage path ####
###

# get taxonomy classifications from NCBI
datax <- read.table("data/kaiju_summary.out", header = T, sep = ";")
rownames(datax) <- datax[,1]
datax <- datax[,-1]
colnames(datax) <- sapply(colnames(datax), function(x){
  cut1 <- strsplit(x, split = "_S[0-9]*_")[[1]][1]
  cut1 <- gsub("output.","",cut1)
  return(cut1)
})

# insert 0 for NA counts
datax[is.na(datax)] <- 0

# make annotations
tax_IDs <- rownames(datax)
taxa_data <- classification(tax_IDs, db = "ncbi")
saveRDS(taxa_data, file = "data/taxa_annotations.rds")
# taxa_data <- readRDS("data/taxa_annotations.rds")

# get annotations for each rank
ranks <- c("superkingdom","phylum","class","order","family","genus","species")
taxa_data_names <- taxid2name(tax_IDs, db = "ncbi")
taxa_data_ranks <- sapply(taxa_data, function(x){
  if(!is.null(dim(x))){
    path <- t(x[match(ranks,x[,2]),"name", drop = FALSE])
    return(path)
  } else {
    return(matrix(rep("NA",length(ranks)), nrow = 1))
  }
}, simplify = TRUE)
taxa_data_ranks <- data.frame(TaxID = colnames(taxa_data_ranks), Name = taxa_data_names, t(taxa_data_ranks))
colnames(taxa_data_ranks) <- c("TaxID", "TaxName", ranks)
taxa_data_ranks[taxa_data_ranks=="NA"] <- NA

# get the lowest level in the lineage path
taxa_ranks <- apply(taxa_data_ranks[,ranks], 1, function(x){
  z <- ranks[!is.na(x)]
  return(ifelse(length(z) > 0, z[length(z)], "no rank"))
}, simplify = TRUE)
taxa_ranks <- unlist(taxa_ranks)
taxa_data_ranks$Rank <- taxa_ranks
saveRDS(taxa_data_ranks, file = "data/taxonomy_ranks_data.rds")