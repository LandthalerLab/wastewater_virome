# Note: Run the entire script :D 

####
# Library ####
####

library(dplyr)
library(umap)
library(ggplot2)
library(Seurat)
library(stringr)

####
# Import Data ####
####

# import embeddings 
orf_embeddings <- read.csv("ProtTrans/orf_embeddings.csv")
orf_embeddings <- orf_embeddings[,-1]
orf_embeddings <- t(orf_embeddings)

# import metadata
metadata <- read.csv("trinity_hmmer_filtered_withBLASTandRNArep.csv")
metadata$Hmm <- sapply(metadata$AllHmm, function(x) {
  strsplit(x, split = " ")[[1]][1]
})

# get embeddings of ORFs analyzed by BLAST/CCtyper
orf_embeddings <- orf_embeddings[metadata$ORF,]

####
# Configure Labels ####
####

# get labels
metadata$HmmClass <- sapply(metadata$Hmm, function(x){
  strsplit(x, split = "_")[[1]][1]
})
temp <- metadata$HmmClass
temp[grepl("Cas|Csc|Csb|Csf|Cse|Csx|Csm|Cmr", temp)] <- str_extract(string = temp[grepl("Cas|Csc|Csb|Csf|Cse|Csx|Csm|Cmr", temp)], 
                                                                    pattern = "(Cas|Csc|Csb|Csf|Cse|Csx|Csm|Cmr)[0-9]+")
temp[metadata$HmmClass=="CasR"] <- "CasR"
metadata$HmmMajorClass <- temp

####
# PCA and UMAP ####
####

# standardize and PCA
orf_embeddings_scaled <- apply(orf_embeddings, 2, scale)
orf_embeddings_pr <- prcomp(orf_embeddings_scaled)
var_perc <- cumsum(orf_embeddings_pr$sdev^2)/sum(orf_embeddings_pr$sdev^2)
plot(1:30, (orf_embeddings_pr$sdev^2)[1:30])
prdim <- 30
orf_embeddings_pr_data <- orf_embeddings_pr$x[,1:prdim]
rownames(orf_embeddings_pr_data) <- rownames(orf_embeddings)

# reduce woth UMAP
set.seed(1)
orf_embeddings_umap <- umap(orf_embeddings_pr_data)
orf_embeddings_umap_data <- data.frame(orf_embeddings_umap$layout)
colnames(orf_embeddings_umap_data) <- c("x", "y")

####
# UMAP ####
####

# get cluster centers
cluster_info_hmmclass <- aggregate(orf_embeddings_umap_data, list(metadata$HmmMajorClass), mean)
colnames(cluster_info_hmmclass) <- c("label", "x", "y")

# plot embeddings
g1 <- ggplot() + 
  geom_point(aes(x = x, y = y, color = as.factor(metadata$HmmMajorClass)), data = orf_embeddings_umap_data) +
  geom_label(data = cluster_info_hmmclass, mapping = aes(x = x, y = y, label = label), segment.color = 'grey50') + 
  guides(fill=guide_legend(ncol=2)) 