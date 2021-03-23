#remove genes that don't have at least 1 cpm in >= 20% of samples in a single project
tic <- as.integer(as.POSIXct(Sys.time()))
library(tidyverse)
args <- commandArgs(TRUE)
#arg 1 is path to directory to filter
#arg 2 is path to file with list of genes to keep

#list files in dir to filter
dir_files <- list.files(args[1], full.names = T)

#create dir for output
dirname <- paste0("./", basename(args[1]), "_gene_filtered")
if(!dir.exists(dirname)) {
  dir.create(dirname)
}

#read in file of genes to keep
genes_to_keep <- read_delim(args[2],
                            delim = "\t",
                            col_names = F) %>% 
  pull(1)

for (filename in dir_files){
  #read file
  file_data <- read_delim(filename, 
                          delim = "\t", 
                          col_names = T, 
                          col_types = cols(.default = "d", gene = "c"))
  #get file basename
  file_basename <- basename(filename)
  #if not .pcl, make .pcl for the sake of sleipnir later on
  file_basename <- gsub("\\.txt", ".pcl", file_basename)
  file_basename <- gsub("\\.tsv", ".pcl", file_basename)
  #remove genes not in genes_to_keep, rewrite file
  file_data <- file_data %>% filter(gene %in% genes_to_keep)
  file_data %>% as.data.frame() %>% 
    write_delim(paste0(dirname, "/", gsub("\\.pcl$", "_cpm_filtered.pcl", file_basename)),
                delim = "\t",
                col_names = T)
}

#script time
toc <- as.integer(as.POSIXct(Sys.time()))
print(paste("The time it took to run this script in minutes was", (toc-tic)/60, sep = " "))
