# library
library(dplyr)

# set arguements
args = commandArgs(trailingOnly=TRUE)

# list of files
list_of_table <- list.files("results/kaiju/wwDNAall_04112022/data/tables/", full.names = TRUE)
list_of_table <- list_of_table[grepl("_kaiju_table.csv",list_of_table)]
# list_of_table <- strsplit(args[1], split = " ")[[1]]

# merge tables
alldata <- NULL
for(i in 1:length(list_of_table)){
  print(list_of_table[i])
  cur_data <- read.table(list_of_table[i], header = T, sep = ";")
  if(i == 1){
    alldata <- cur_data
  } else {
    alldata <- alldata %>% full_join(cur_data)
  }
}

