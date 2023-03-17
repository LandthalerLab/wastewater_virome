# library
library(dplyr)
library(taxizedb)
library(taxonomizr)
library(seqinr)
library(rentrez)
library(XML)
library(xlsx)
library(ggplot2)
# options(java.parameters = "-Xmx100000m")

####
# Import Clustering results of CD-HIT ####
####

clstr <- read.csv("results/cctyper/RNADNA_combined/findNonredundant_ORF/results.clstr", sep = "\t", row.names = NULL, header = FALSE, stringsAsFactors = FALSE)

####
# Change data ####
####

clstr2 <- clstr
n = nrow(clstr)
x = 0
numbers_only <- function(x) !grepl("\\D", x)
for (row in c(1:n)) {
  if (numbers_only(clstr2[row,1]) == TRUE) {
    clstr2[row,1] <- x}
  else {NULL}
  x <- clstr2[row,1]
}
clstr.sums <- data.frame(dplyr::count(clstr2,V1))
clstr.sums <- clstr.sums[order(clstr.sums$n, decreasing = TRUE),]
clstr2$V1 <- gsub(">", "", clstr2$V1)
clstr2$V2 <- lapply(clstr2$V2, function(x) {
  temp <- strsplit(x, split = ">")[[1]][2]
  temp <- strsplit(temp, split = "\\.")[[1]][1]
  temp
})
clstr2 <- clstr2[!is.na(clstr2$V2),]
# clstr2$V2 <- gsub("^","TRINITY_DN",clstr2$V2)
colnames(clstr2) <- c("Cluster", "ORF")

####
# Change data ####
####

write.table(as.matrix(clstr2), file = "results/cctyper/RNADNA_combined/findNonredundant_ORF/clusters_ORF.txt", sep = "\t")


