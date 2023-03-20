library(rBLAST)

# get arguements
args <- commandArgs(trailingOnly=TRUE)

# get the reads and filename
file <- args[1]
output_file <- paste0("results/", basename(file), ".rds")

# run blast
bl <- blast(db = "/fast/AG_Landthaler/genomes/BLAST/nr/fasta/nr.fasta", type = "blastp")
seq <- readAAStringSet(file)
cl <- predict(bl, seq)

# save file
saveRDS(cl, file = output_file)
