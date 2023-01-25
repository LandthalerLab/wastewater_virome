# library
library(dplyr)

# set arguements
args = commandArgs(trailingOnly=TRUE)

# list of files
# list_of_table <- list.files(args[1], full.names = TRUE)
# list_of_table <- list_of_table[grepl("_kaiju_table.out",list_of_table)]
list_of_table <- strsplit(args[1], split = " ")[[1]]

# merge tables
alldata <- NULL
for(i in 1:length(list_of_table)){
  print(list_of_table[i])
  cur_data <- read.table(list_of_table[i], header = T, sep = ";")
  # cur_data$taxon_id[cur_data$taxon_name == "unclassified"] <- "unclassified"
  # cur_data$taxon_id[grepl("cannot be assigned", cur_data$taxon_name)] <- "Multimapped"
  # cur_data$taxon_id <- paste0(cur_data$taxon_id, "_", cur_data$taxon_name) 
  # cur_data <- cur_data[,c("taxon_id","reads")]
  # colnames(cur_data) <- c("ta", list_of_table[i])
  if(i == 1){
    alldata <- cur_data
  } else {
    alldata <- alldata %>% full_join(cur_data)
  }
}

write.table(alldata, args[2], quote = FALSE, row.names = FALSE, sep = ";")
