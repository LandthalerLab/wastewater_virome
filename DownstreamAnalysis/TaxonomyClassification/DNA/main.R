# library
library(dplyr)
library(ComplexHeatmap)
library(ggplot2)
library(reshape)
library(taxizedb)
library(debrowser)
library(heatmaply)
library(circlize)
library(tidyverse)
library(tibble)
library(ggfortify)
library(ggrepel)

###
# Import Data ####
###

# import kaiju 
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

# import dedup stats
dedup_data <- read.table("data/all_dedup_summary.out", header = T, sep = "\t", quote = "")

# import taxonomy data
taxa_data_ranks <- readRDS("data/taxonomy_ranks_data.rds")

###
# Low Count Filtering ####
###

# quality filtering and dont filter out some features
datax_max_count <- apply(datax, 1, max)
non_filtered_tax <- c("272636","10524")
datax_highcount <- datax[datax_max_count > 10 | rownames(datax) %in% non_filtered_tax, ]

###
# Process Data ####
###

# normalize count matrix, hellinger norm without unassigned
datax_norm <- apply(datax_highcount[!rownames(datax_highcount) %in% c("0","1","2","10239","2759","2157"),], 2, function(x){
  return(sqrt(((x)/sum(x))*1000000))
})

###
# PCA ####
###

datax_de <- datax_norm
datax_de <- apply(datax_de,1,scale)
datax_pr <- datax_de
prtemp <- prcomp(datax_pr)
datax_pr_x <- prtemp$x
rownames(datax_pr_x) <- colnames(datax_norm)
var_perc <- cumsum(prtemp$sdev^2)/sum(prtemp$sdev^2)
ggplot(data.frame(datax_pr_x), aes(x = PC1, y = PC2, label = rownames(datax_pr_x))) + geom_text()
ggplot_pca <- ggplot2::autoplot(prtemp, data = datax_pr) +
  theme_classic() + 
  theme(axis.ticks.x = element_blank(), axis.text.x = element_blank(), 
        axis.ticks.y = element_blank(), axis.text.y = element_blank())
ggplot_pca

###
# Heatmap ####
###

###
## Genus-wise counts for days ####
###

# data
datax_assigned <- datax_highcount[!rownames(datax_highcount) %in% c("0","1","2","10239","2759","2157"),]

# aggregate code from emanuel
mapping_table <- as_tibble(datax_assigned)
mapping_table <- as.data.frame(t(mapping_table))
mapping_table$day <- gsub(pattern = "^DNA_([0-9]{6})_.*", replacement = "\\1", rownames(mapping_table), perl = TRUE)
mapping_table <- mapping_table %>% group_by(day) %>% summarize_all(sum) %>% column_to_rownames(var="day")
mapping_table <- as.data.frame(t(mapping_table))
rownames(mapping_table) <- rownames(datax_assigned)

# size normalized count data
mapping_table <- apply(mapping_table, 2, function(x){
  return((x/sum(x))*1000000)
})

# aggregate by family, separate for superkingdom
taxa_data_ranks_subset <- taxa_data_ranks[rownames(mapping_table),]
genus <- taxa_data_ranks_subset$genus
genus[is.na(genus)] <- "None"
norm_datax_orders <- data.frame(mapping_table, genus = genus) %>%
  group_by(genus) %>% summarise_at(1:ncol(mapping_table), sum, na.rm = TRUE) 
taxa_data_ranks_subset$genus[is.na(taxa_data_ranks_subset$genus)] <- "None"
taxa_data_ranks_subset_group <- data.frame(taxa_data_ranks_subset) %>%
  group_by(genus) %>% summarise(superkingdom = superkingdom[1])
norm_datax_orders <- as.data.frame(norm_datax_orders)
norm_datax_orders <- norm_datax_orders[!taxa_data_ranks_subset_group$genus == "None",]
taxa_data_ranks_subset_group <- taxa_data_ranks_subset_group[!taxa_data_ranks_subset_group$genus == "None",]

# fix row names
colnames(norm_datax_orders) <- c("genus",colnames(mapping_table))

# order each superkingdom by mean count abundance
# choose top 20 family from each kingdom
norm_datax_orders$Total <- rowSums(norm_datax_orders[,-1])
norm_datax_orders$superkingdom <- taxa_data_ranks_subset_group$superkingdom
norm_datax_orders_top <- as_tibble(norm_datax_orders) %>% 
  group_by(superkingdom) %>% 
  slice_max(order_by = Total, n = 25)
taxa_data_ranks_subset_group_top <- taxa_data_ranks_subset_group[match(norm_datax_orders_top$genus,taxa_data_ranks_subset_group$genus), ]

# fix row names
norm_datax_orders_top <- as.data.frame(norm_datax_orders_top)
rownames(norm_datax_orders_top) <- norm_datax_orders_top[,1]
norm_datax_orders_top <- norm_datax_orders_top[,-1]
norm_datax_orders_top <- norm_datax_orders_top[,!colnames(norm_datax_orders_top) %in% c("Total","superkingdom")]

# log normalize 
norm_datax_orders_top_log <- log10(norm_datax_orders_top + 1)

# bacteria
heatmap_data <- norm_datax_orders_top_log[taxa_data_ranks_subset_group_top$superkingdom=="Bacteria",]
col_fun = colorRamp2(c(3.5, 4, 4.8), c("dodgerblue", "#fbfcbd", "red"))
bacteria_heatmap <- Heatmap(as.matrix(heatmap_data), col = col_fun,
                            cluster_columns = FALSE, cluster_rows = TRUE, show_row_dend = FALSE, show_column_dend = FALSE, 
                            show_row_names = TRUE, show_column_names = TRUE, column_names_rot = 45, column_title = "", 
                            # row_split = taxa_data_ranks_subset_group_top$superkingdom,
                            row_title_rot = 0, row_names_side = "left", row_names_gp = grid::gpar(fontsize = 9), 
                            column_names_gp = grid::gpar(fontsize = 15), heatmap_legend_param = list(title = "Log10.Count.\n(Bacteria)"))

# eukaryota
heatmap_data <- norm_datax_orders_top_log[taxa_data_ranks_subset_group_top$superkingdom=="Eukaryota",]
col_fun = colorRamp2(c(1, 3, 4.75), c("dodgerblue", "#fbfcbd", "red"))
eukaryota_heatmap <- Heatmap(as.matrix(heatmap_data), col = col_fun,
                             cluster_columns = FALSE, cluster_rows = TRUE, show_row_dend = FALSE, show_column_dend = FALSE, 
                             show_row_names = TRUE, show_column_names = TRUE, column_names_rot = 45, column_title = "", 
                             # row_split = taxa_data_ranks_subset_group_top$superkingdom,
                             row_title_rot = 0, row_names_side = "left", row_names_gp = grid::gpar(fontsize = 9), 
                             column_names_gp = grid::gpar(fontsize = 15), heatmap_legend_param = list(title = "Log10.Count. \n (Eukaryota)"))

# viruses
heatmap_data <- norm_datax_orders_top_log[taxa_data_ranks_subset_group_top$superkingdom=="Viruses",]
col_fun = colorRamp2(c(1, 3, 5), c("dodgerblue", "#fbfcbd", "red"))
viruses_heatmap <- Heatmap(as.matrix(heatmap_data), col = col_fun,
                           cluster_columns = FALSE, cluster_rows = TRUE, show_row_dend = FALSE, show_column_dend = FALSE, 
                           show_row_names = TRUE, show_column_names = TRUE, column_names_rot = 45, column_title = "", 
                           # row_split = taxa_data_ranks_subset_group_top$superkingdom,
                           row_title_rot = 0, row_names_side = "left", row_names_gp = grid::gpar(fontsize = 9), 
                           column_names_gp = grid::gpar(fontsize = 15), heatmap_legend_param = list(title = "Log10.Count. \n (Viruses)"))

# archaea
heatmap_data <- norm_datax_orders_top_log[taxa_data_ranks_subset_group_top$superkingdom=="Archaea",]
col_fun = colorRamp2(c(1, 2, 4), c("dodgerblue", "#fbfcbd", "red"))
archaea_heatmap <- Heatmap(as.matrix(heatmap_data), col = col_fun,
                           cluster_columns = FALSE, cluster_rows = TRUE, show_row_dend = FALSE, show_column_dend = FALSE, 
                           show_row_names = TRUE, show_column_names = TRUE, column_names_rot = 45, column_title = "", 
                           # row_split = taxa_data_ranks_subset_group_top$superkingdom,
                           row_title_rot = 0, row_names_side = "left", row_names_gp = grid::gpar(fontsize = 9), 
                           column_names_gp = grid::gpar(fontsize = 15), heatmap_legend_param = list(title = "Log10.Count. \n (Archaea)"))

# get selected taxonomies
selected <- c("272636","10524")
selected_intersect <- intersect(selected, rownames(mapping_table))
mapping_table_zoom <- mapping_table[selected_intersect,]

# 0-1 normalize 
mapping_table_zoom_norm <- t(apply(mapping_table_zoom, 1, function(x) {
  (x -min(x))/(max(x)- min(x))
}))

# visualize
rownames(mapping_table_zoom_norm) <- taxa_data_ranks$TaxName[match(rownames(mapping_table_zoom), taxa_data_ranks$TaxID)]
col_fun = colorRamp2(c(0, 0.5, 1), c("dodgerblue", "#fbfcbd", "red"))
selectedHeatmap <- Heatmap(as.matrix(mapping_table_zoom_norm), col = col_fun,
                           cluster_columns = FALSE, cluster_rows = FALSE, show_row_dend = FALSE, show_column_dend = FALSE, 
                           show_row_names = TRUE, show_column_names = TRUE, column_names_rot = 45, column_title = "", 
                           row_names_rot = 0, row_names_side = "left", row_names_gp = grid::gpar(fontsize = 9), 
                           column_names_gp = grid::gpar(fontsize = 12),  heatmap_legend_param = list(title = "Norm. Count"),
                           # show_heatmap_legend = FALSE
                           height = nrow(mapping_table_zoom_norm)*unit(3, "mm"))

# Return Heatmap
bacteria_heatmap %v% eukaryota_heatmap %v% viruses_heatmap %v% archaea_heatmap %v% selectedHeatmap
