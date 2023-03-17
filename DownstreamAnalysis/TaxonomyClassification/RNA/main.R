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
# Import data ####
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

# remove feb 2022 and july 2021
datax <- datax[,!grepl("2107|2202", colnames(datax))]

# insert 0 for NA counts
datax[is.na(datax)] <- 0

# import dedup stats
dedup_data <- read.table("data/all_dedup_summary.out", header = T, sep = "\t", quote = "")
dedup_data <- dedup_data[!grepl("2107|2202", dedup_data$Sample),]

# import taxonomy data
taxa_data_ranks <- readRDS("data/taxonomy_ranks_data.rds")

###
# Statistics on Duplication, Assignment and SuperKingdoms ####
###

# get duplicated read counts
dedup_rate <- dedup_data$Total.Reads - sapply(dedup_data$NoDup.Reads, function(x) return(as.numeric(strsplit(x, split = " ")[[1]][1])))

# assigned, unassigned and duplicate read counts
superkingdom <- taxa_data_ranks$superkingdom
datax_stats <- data.frame(sample = colnames(datax),
                          sizeFactor = colSums(datax),
                          Unassigned = colSums(datax[rownames(datax) == "0" | rownames(datax) == "1",]),
                          Viruses = colSums(datax[superkingdom=="Viruses",], na.rm = TRUE),
                          Bacteria = colSums(datax[superkingdom=="Bacteria",], na.rm = TRUE),
                          Archaea = colSums(datax[superkingdom=="Archaea",], na.rm = TRUE),
                          Eukaryota = colSums(datax[superkingdom=="Eukaryota",], na.rm = TRUE),
                          Duplicate = dedup_rate)
datax_stats$Assigned <- rowSums(datax_stats[,c("Viruses","Bacteria","Archaea","Eukaryota")])
datax_stats$NonDuplicate <- rowSums(datax_stats[,c("Viruses","Bacteria","Archaea","Eukaryota","Unassigned")])

# duplicate
assignment_stats <- melt(datax_stats, measure.vars = c("NonDuplicate", "Duplicate"), variable_name = "Status")
ggplot(data=assignment_stats, aes(x=gsub("WW_","", sample), y=value, fill=Status)) +
  geom_bar(stat="identity") + xlab("Samples") + ylab("# of Reads") + 
  guides(fill=guide_legend(title="Read Status")) + 
  theme_classic() + 
  theme(axis.ticks.x = element_blank(), 
        axis.ticks.y = element_blank(),
        axis.text.x=element_text(angle=45, hjust=1, vjust = 1, size = rel(0.7))) + 
  scale_y_continuous(breaks = c(0,10000000, 20000000, 30000000), limits = c(0,30010000), labels = scales::comma)  

# superkingdom and unassigned
kingdom_stats <- datax_stats[,c("sample", "Viruses", "Bacteria", "Archaea", "Eukaryota", "Unassigned")]
kingdom_stats <- melt(kingdom_stats, measure.vars = c("Viruses", "Bacteria", "Archaea", "Eukaryota", "Unassigned"), variable_name = "SuperKingdom")
kingdom_stats <- as_tibble(kingdom_stats) %>% group_by(sample) %>% mutate(totalreads = sum(value)) %>% mutate(prop = value/totalreads)
ggplot(data=kingdom_stats, aes(x=gsub("WW_","", sample),
                               y = prop,
                               fill=factor(SuperKingdom, levels = c("Unassigned", "Bacteria", "Eukaryota", "Viruses", "Archaea")))) +
  geom_bar(stat="identity") + xlab("Samples") + ylab("# of Reads") + 
  theme_classic() + 
  theme(axis.ticks.x = element_blank(), 
        axis.ticks.y = element_blank(),
        axis.text.x=element_text(angle=45, hjust=1, vjust = 1, size = rel(0.7))) + 
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) + 
  guides(fill=guide_legend(title="Super Kingdom")) 

###
# Low Filter Count ####
###

# quality filtering
datax_max_count <- apply(datax, 1, max)
datax_highcount <- datax[datax_max_count > 10, ]

# normalize count matrix, hellinger norm without unassigned
datax_norm <- apply(datax_highcount[!rownames(datax_highcount) %in% c("0","1","2","10239","2759","2157"),], 2, function(x){
  return(sqrt(((x)/sum(x))*1000000))
})

###
# Outlier Detection ####
###

###
## Outlying Annotations ####
###

top_n <- 3
outlier_annot_list <- NULL
top_annot <- head(rownames(prtemp$rotation[order(prtemp$rotation[,2], decreasing = TRUE),]),top_n)
bottom_annot <- head(rownames(prtemp$rotation[order(prtemp$rotation[,1], decreasing = TRUE),]),top_n)
outlier_annot_list <- c(outlier_annot_list, c(top_annot, bottom_annot))

###
## PCA ####
###

ggplot_pca <- ggplot2::autoplot(prtemp, data = datax_pr, loadings = TRUE, loadings.label=TRUE, 
                                loadings.colour = "red", loadings.label.repel=T, loadings.alpha = 0.2) + 
  theme_classic() + 
  theme(axis.ticks.x = element_blank(), axis.text.x = element_blank(), 
        axis.ticks.y = element_blank(), axis.text.y = element_blank())
ggplot_pca$layers[[2]]$data <- ggplot_pca$layers[[2]]$data[outlier_annot_list,]
ggplot_pca$layers[[3]]$data <- ggplot_pca$layers[[3]]$data[outlier_annot_list,]
ggplot_pca

## Outlying Annotations vs Outliers ####

# quality filtering
datax_mean_count <- apply(datax, 1, mean)
datax_highcount2 <- datax[datax_mean_count > 10, ]

# normalize count matrix, hellinger norm without unassigned
datax_norm <- apply(datax_highcount2[!rownames(datax_highcount2) %in% c("0","1","2","10239","2759","2157"),], 2, function(x){
  return(sqrt(((x)/sum(x))*1000000))
})

### get annotations different in new filtering and old filtering
datax_highcount_new <- datax_highcount[!rownames(datax_highcount) %in% c("0","1","2","10239","2759","2157"),]
datax_highcount2_new <- datax_highcount2[!rownames(datax_highcount2) %in% c("0","1","2","10239","2759","2157"),]

removed_annotations <- setdiff(rownames(datax_highcount_new), rownames(datax_highcount2_new))
removed_annotations_taxID <- ifelse(rownames(datax_highcount_new) %in% removed_annotations,
                                    "Removed","Preserved")

datax_removedvspreserved <- aggregate(datax_highcount_new, list(removed_annotations_taxID), sum)
datax_removedvspreserved <- melt(datax_removedvspreserved)
datax_removedvspreserved$variable <- gsub("WW_","",datax_removedvspreserved$variable)
datax_removed <- datax_removedvspreserved[datax_removedvspreserved$Group.1=="Removed",]
top_outliers <- head(as.character(datax_removed$variable)[order(datax_removed$value, decreasing = TRUE)], 10)
datax_removed$top_outliers <- ifelse(datax_removed$variable %in% top_outliers, "Outlier", "Non-Outlier") 
ggplot(data=datax_removed, 
       aes(x = reorder(variable, -value), y = value, fill = top_outliers)) +
  geom_bar(stat="identity") + 
  theme(axis.text.x=element_text(angle=45, hjust=1, vjust = 1, size = rel(0.7))) + 
  guides(fill=guide_legend(title="")) + 
  xlab("Samples") + ylab("# of Reads") +
  scale_y_continuous(breaks = c(0,50000, 100000, 150000), limits = c(0,150010), labels = scales::comma)  

###
# Heatmap ####
###

###
## Genus-wise counts for Months ####
###

# data
datax_assigned <- datax_highcount[!rownames(datax_highcount) %in% c("0","1","2","10239","2759","2157"),]

# aggregate code from emanuel
mapping_table <- as_tibble(datax_assigned)
mapping_table <- mapping_table %>% select(-contains(c("mix_1", "mix_2", "BIOM", "PCRinh", "captu", "cer_")))
mapping_table <- as.data.frame(t(mapping_table))
mapping_table$month <- gsub(pattern = "^WW_([0-9]{4}).*", replacement = "\\1", rownames(mapping_table), perl = TRUE)
mapping_table <- mapping_table %>% group_by(month) %>% summarize_all(sum) %>% column_to_rownames(var="month")
mapping_table <- as.data.frame(t(mapping_table))
rownames(mapping_table) <- rownames(datax_assigned)
mapping_table <- mapping_table[,!colnames(mapping_table) %in% c("2107","2202")]

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
                            row_title_rot = 0, row_names_side = "left", row_names_gp = grid::gpar(fontsize = 9), 
                            column_names_gp = grid::gpar(fontsize = 15), heatmap_legend_param = list(title = "Log10.Count.\n(Bacteria)"))

# eukaryota
heatmap_data <- norm_datax_orders_top_log[taxa_data_ranks_subset_group_top$superkingdom=="Eukaryota",]
col_fun = colorRamp2(c(1, 3, 4.75), c("dodgerblue", "#fbfcbd", "red"))
eukaryota_heatmap <- Heatmap(as.matrix(heatmap_data), col = col_fun,
                             cluster_columns = FALSE, cluster_rows = TRUE, show_row_dend = FALSE, show_column_dend = FALSE, 
                             show_row_names = TRUE, show_column_names = TRUE, column_names_rot = 45, column_title = "", 
                             row_title_rot = 0, row_names_side = "left", row_names_gp = grid::gpar(fontsize = 9), 
                             column_names_gp = grid::gpar(fontsize = 15), heatmap_legend_param = list(title = "Log10.Count. \n (Eukaryota)"))

# viruses
heatmap_data <- norm_datax_orders_top_log[taxa_data_ranks_subset_group_top$superkingdom=="Viruses",]
col_fun = colorRamp2(c(1, 3, 5), c("dodgerblue", "#fbfcbd", "red"))
viruses_heatmap <- Heatmap(as.matrix(heatmap_data), col = col_fun,
                           cluster_columns = FALSE, cluster_rows = TRUE, show_row_dend = FALSE, show_column_dend = FALSE, 
                           show_row_names = TRUE, show_column_names = TRUE, column_names_rot = 45, column_title = "", 
                           row_title_rot = 0, row_names_side = "left", row_names_gp = grid::gpar(fontsize = 9), 
                           column_names_gp = grid::gpar(fontsize = 15), heatmap_legend_param = list(title = "Log10.Count. \n (Viruses)"))

# archaea
heatmap_data <- norm_datax_orders_top_log[taxa_data_ranks_subset_group_top$superkingdom=="Archaea",]
col_fun = colorRamp2(c(1, 2, 4), c("dodgerblue", "#fbfcbd", "red"))
archaea_heatmap <- Heatmap(as.matrix(heatmap_data), col = col_fun,
                           cluster_columns = FALSE, cluster_rows = TRUE, show_row_dend = FALSE, show_column_dend = FALSE, 
                           show_row_names = TRUE, show_column_names = TRUE, column_names_rot = 45, column_title = "", 
                           row_title_rot = 0, row_names_side = "left", row_names_gp = grid::gpar(fontsize = 9), 
                           column_names_gp = grid::gpar(fontsize = 15), heatmap_legend_param = list(title = "Log10.Count. \n (Archaea)"))

# get selected taxonomies
selected <- c("1758882","103722","10941","1803956","146500","1239565","12239")
selected_intersect <- intersect(selected, rownames(mapping_table))
mapping_table_zoom <- mapping_table[selected_intersect,]
mapping_table_zoom <- mapping_table_zoom[c("12239", "103722", "1239565", "1758882"),]

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
                           height = nrow(mapping_table_zoom)*unit(3, "mm"))

# return heatmap
bacteria_heatmap %v% eukaryota_heatmap %v% viruses_heatmap %v% archaea_heatmap %v% selectedHeatmap

###
## Genus-wise counts for Samples ####
###

# data
mapping_table <- datax_highcount[!rownames(datax_highcount) %in% c("0","1","2","10239","2759","2157"),]
colnames(mapping_table) <- gsub("WW_","", colnames(mapping_table))

# aggregate code from emanuel
mapping_table <- mapping_table[,!grepl("2107|2202",colnames(mapping_table))]

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
                             column_names_gp = grid::gpar(fontsize = 15), heatmap_legend_param = list(title = "Log10.Count.\n(Eukaryota)"))

# viruses
heatmap_data <- norm_datax_orders_top_log[taxa_data_ranks_subset_group_top$superkingdom=="Viruses",]
col_fun = colorRamp2(c(1, 3, 5), c("dodgerblue", "#fbfcbd", "red"))
viruses_heatmap <- Heatmap(as.matrix(heatmap_data), col = col_fun,
                           cluster_columns = FALSE, cluster_rows = TRUE, show_row_dend = FALSE, show_column_dend = FALSE, 
                           show_row_names = TRUE, show_column_names = TRUE, column_names_rot = 45, column_title = "", 
                           row_title_rot = 0, row_names_side = "left", row_names_gp = grid::gpar(fontsize = 9), 
                           column_names_gp = grid::gpar(fontsize = 15), heatmap_legend_param = list(title = "Log10.Count.\n(Viruses)"))

# archaea
heatmap_data <- norm_datax_orders_top_log[taxa_data_ranks_subset_group_top$superkingdom=="Archaea",]
col_fun = colorRamp2(c(1, 2, 4), c("dodgerblue", "#fbfcbd", "red"))
archaea_heatmap <- Heatmap(as.matrix(heatmap_data), col = col_fun,
                           cluster_columns = FALSE, cluster_rows = TRUE, show_row_dend = FALSE, show_column_dend = FALSE, 
                           show_row_names = TRUE, show_column_names = TRUE, column_names_rot = 45, column_title = "", 
                           row_title_rot = 0, row_names_side = "left", row_names_gp = grid::gpar(fontsize = 9), 
                           column_names_gp = grid::gpar(fontsize = 15), heatmap_legend_param = list(title = "Log10.Count.\n(Archaea)"))

# get selected taxonomies
selected <- c("1758882","103722","10941","1803956","146500","1239565","12239")
selected_intersect <- intersect(selected, rownames(mapping_table))
mapping_table_zoom <- mapping_table[selected_intersect,]
mapping_table_zoom <- mapping_table_zoom[c("12239", "103722", "1239565", "1758882"),]

# 0-1 normalize 
# mapping_table_zoom_log <- log10(mapping_table_zoom + 1)
mapping_table_zoom_norm <- t(apply(mapping_table_zoom, 1, function(x) {
  (x -min(x))/(max(x)- min(x))
}))

# visualize
rownames(mapping_table_zoom_norm) <- taxa_data_ranks$TaxName[match(rownames(mapping_table_zoom), taxa_data_ranks$TaxID)]
col_fun = colorRamp2(c(0, 0.5, 1), c("dodgerblue", "#fbfcbd", "red"))
selectedHeatmap <- Heatmap(as.matrix(mapping_table_zoom_norm), col = col_fun,
                           cluster_columns = FALSE, cluster_rows = FALSE, show_row_dend = FALSE, show_column_dend = FALSE, 
                           show_row_names = TRUE, show_column_names = TRUE, column_names_rot = 45, column_title = "", 
                           row_names_rot = 0, row_names_side = "left",
                           column_names_gp = grid::gpar(fontsize = 5),  row_names_gp = grid::gpar(fontsize = 9),  
                           heatmap_legend_param = list(title = "Norm. Count"),
)

# return heatmap
bacteria_heatmap %v% eukaryota_heatmap %v% viruses_heatmap %v% archaea_heatmap %v% selectedHeatmap

###
## Family-wise counts for Months ####
###

# data
datax_assigned <- datax_highcount[!rownames(datax_highcount) %in% c("0","1","2","10239","2759","2157"),]

# aggregate code from emanuel
mapping_table <- as_tibble(datax_assigned)
mapping_table <- mapping_table %>% select(-contains(c("mix_1", "mix_2", "BIOM", "PCRinh", "captu", "cer_")))
mapping_table <- as.data.frame(t(mapping_table))
mapping_table$month <- gsub(pattern = "^WW_([0-9]{4}).*", replacement = "\\1", rownames(mapping_table), perl = TRUE)
mapping_table <- mapping_table %>% group_by(month) %>% summarize_all(sum) %>% column_to_rownames(var="month")
mapping_table <- as.data.frame(t(mapping_table))
rownames(mapping_table) <- rownames(datax_assigned)
mapping_table <- mapping_table[,!colnames(mapping_table) %in% c("2107","2202")]

# size normalized count data
mapping_table <- apply(mapping_table, 2, function(x){
  return((x/sum(x))*1000000)
})

# aggregate by family, separate for superkingdom
taxa_data_ranks_subset <- taxa_data_ranks[rownames(mapping_table),]
family <- taxa_data_ranks_subset$family
family[is.na(family)] <- "None"
norm_datax_orders <- data.frame(mapping_table, family = family) %>%
  group_by(family) %>% summarise_at(1:ncol(mapping_table), sum, na.rm = TRUE) 
taxa_data_ranks_subset$family[is.na(taxa_data_ranks_subset$family)] <- "None"
taxa_data_ranks_subset_group <- data.frame(taxa_data_ranks_subset) %>%
  group_by(family) %>% summarise(superkingdom = superkingdom[1])
norm_datax_orders <- as.data.frame(norm_datax_orders)
norm_datax_orders <- norm_datax_orders[!taxa_data_ranks_subset_group$family == "None",]
taxa_data_ranks_subset_group <- taxa_data_ranks_subset_group[!taxa_data_ranks_subset_group$family == "None",]

# fix row names
colnames(norm_datax_orders) <- c("family",colnames(mapping_table))

# order each superkingdom by mean count abundance
# choose top 20 family from each kingdom
norm_datax_orders$Total <- rowSums(norm_datax_orders[,-1])
norm_datax_orders$superkingdom <- taxa_data_ranks_subset_group$superkingdom
norm_datax_orders_top <- as_tibble(norm_datax_orders) %>% 
  group_by(superkingdom) %>% 
  slice_max(order_by = Total, n = 25)
taxa_data_ranks_subset_group_top <- taxa_data_ranks_subset_group[match(norm_datax_orders_top$family,taxa_data_ranks_subset_group$family), ]

# fix row names
norm_datax_orders_top <- as.data.frame(norm_datax_orders_top)
rownames(norm_datax_orders_top) <- norm_datax_orders_top[,1]
norm_datax_orders_top <- norm_datax_orders_top[,-1]
norm_datax_orders_top <- norm_datax_orders_top[,!colnames(norm_datax_orders_top) %in% c("Total","superkingdom")]

# log normalize 
norm_datax_orders_top_log <- log10(norm_datax_orders_top + 1)

# log transform and visualize for each superkingdom then merge

# bacteria
heatmap_data <- norm_datax_orders_top_log[taxa_data_ranks_subset_group_top$superkingdom=="Bacteria",]
col_fun = colorRamp2(c(3.5, 4, 4.8), c("dodgerblue", "#fbfcbd", "red"))
bacteria_heatmap <- Heatmap(as.matrix(heatmap_data), col = col_fun,
                            cluster_columns = FALSE, cluster_rows = TRUE, show_row_dend = FALSE, show_column_dend = FALSE, 
                            show_row_names = TRUE, show_column_names = TRUE, column_names_rot = 45, column_title = "", 
                            row_title_rot = 0, row_names_side = "left", row_names_gp = grid::gpar(fontsize = 9), 
                            column_names_gp = grid::gpar(fontsize = 15), heatmap_legend_param = list(title = "Log10.Count.\n(Bacteria)"))

# eukaryota
heatmap_data <- norm_datax_orders_top_log[taxa_data_ranks_subset_group_top$superkingdom=="Eukaryota",]
col_fun = colorRamp2(c(1, 2.2, 4.10), c("dodgerblue", "#fbfcbd", "red"))
eukaryota_heatmap <- Heatmap(as.matrix(heatmap_data), col = col_fun,
                             cluster_columns = FALSE, cluster_rows = TRUE, show_row_dend = FALSE, show_column_dend = FALSE, 
                             show_row_names = TRUE, show_column_names = TRUE, column_names_rot = 45, column_title = "", 
                             row_title_rot = 0, row_names_side = "left", row_names_gp = grid::gpar(fontsize = 9), 
                             column_names_gp = grid::gpar(fontsize = 15), heatmap_legend_param = list(title = "Log10.Count. \n (Eukaryota)"))

# viruses
heatmap_data <- norm_datax_orders_top_log[taxa_data_ranks_subset_group_top$superkingdom=="Viruses",]
col_fun = colorRamp2(c(1, 3, 5.20), c("dodgerblue", "#fbfcbd", "red"))
viruses_heatmap <- Heatmap(as.matrix(heatmap_data), col = col_fun,
                           cluster_columns = FALSE, cluster_rows = TRUE, show_row_dend = FALSE, show_column_dend = FALSE, 
                           show_row_names = TRUE, show_column_names = TRUE, column_names_rot = 45, column_title = "", 
                           row_title_rot = 0, row_names_side = "left", row_names_gp = grid::gpar(fontsize = 9), 
                           column_names_gp = grid::gpar(fontsize = 15), heatmap_legend_param = list(title = "Log10.Count. \n (Viruses)"))

# archaea
heatmap_data <- norm_datax_orders_top_log[taxa_data_ranks_subset_group_top$superkingdom=="Archaea",]
col_fun = colorRamp2(c(0, 1, 3.9), c("dodgerblue", "#fbfcbd", "red"))
archaea_heatmap <- Heatmap(as.matrix(heatmap_data), col = col_fun,
                           cluster_columns = FALSE, cluster_rows = TRUE, show_row_dend = FALSE, show_column_dend = FALSE, 
                           show_row_names = TRUE, show_column_names = TRUE, column_names_rot = 45, column_title = "", 
                           row_title_rot = 0, row_names_side = "left", row_names_gp = grid::gpar(fontsize = 9), 
                           column_names_gp = grid::gpar(fontsize = 15), heatmap_legend_param = list(title = "Log10.Count. \n (Archaea)"))

# get selected taxonomies
selected <- c("1758882","103722","10941","1803956","146500","1239565","12239")
selected_intersect <- intersect(selected, rownames(mapping_table))
mapping_table_zoom <- mapping_table[selected_intersect,]
mapping_table_zoom <- mapping_table_zoom[c("12239", "103722", "1239565", "1758882"),]

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
                           height = nrow(mapping_table_zoom_norm)*unit(3, "mm"))

# return heatmap
bacteria_heatmap %v% eukaryota_heatmap %v% viruses_heatmap %v% archaea_heatmap %v% selectedHeatmap

