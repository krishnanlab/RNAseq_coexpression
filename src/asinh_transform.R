#this script transforms all data values with the 
#inverse hyperbolic sine function (asinh)
tic <- as.integer(as.POSIXct(Sys.time()))
library(tidyverse)

args <- commandArgs(TRUE)
#args[1] is directory with files
directory <- args[1]
#list files from dir to transform
file_list <- list.files(directory, full.names = T)

dirname <- paste0("./", basename(directory), "_asinh_transformed")
if(!dir.exists(dirname)) {
  dir.create(dirname)
}

for (file_name in file_list){
  gene_data <- read_delim(file_name, 
                          delim = "\t", 
                          col_names = T, 
                          col_types = cols(.default = "d", gene = "c")) %>% 
    column_to_rownames("gene") %>% 
    as.data.frame()
  #get basename for writing file
  base_filename <- basename(file_name)
  asinh(gene_data) %>% 
    rownames_to_column("gene") %>%
    write_delim(paste0(dirname, "/", gsub("\\.pcl$", "_gene-type_filtered.pcl", file_basename)), 
                delim = "\t", col_names = T)
}

#script time
toc <- as.integer(as.POSIXct(Sys.time()))
print(paste("The time it took to run this script in minutes was", (toc-tic)/60, sep = " "))

