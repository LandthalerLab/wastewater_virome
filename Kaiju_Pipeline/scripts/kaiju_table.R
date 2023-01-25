# library
library(dplyr)

# set arguements
args = commandArgs(trailingOnly=TRUE)

# import kaiju file
data_file <- args[1]
# datax <- read.table(data_file, header=F, sep='\t', fill=T, col.names=paste0('V', 1:8))
datax <- read.delim2(data_file, header=F, fill=T, col.names=paste0('V', 1:7))

# get accessions of reads
# single_index <- paste0(datax[,3],",") ==datax[,5]
# datax_single <- datax[single_index, ]
# datax_single <- cbind(datax_single, 
#                      sapply(datax_single[,6], function(x){
#                        strsplit(x, split = ",")[[1]][1]
#                      }))
# datax_single[,3] <- paste0(datax_single[,3], "_", gsub(",","", datax_single[,9]))
# datax[single_index, ] <- datax_single

# attach lineage path
# datax[,8] <- gsub(" ","_",datax[,8])
# datax[,3] <- paste0(datax[,3], "_", datax[,8])

# make table
datax_table <- as.data.frame(table(datax[,3]))
colnames(datax_table) <- c("TaxID", args[1])
write.table(datax_table, args[2], quote = FALSE, row.names = FALSE, sep = ";")
