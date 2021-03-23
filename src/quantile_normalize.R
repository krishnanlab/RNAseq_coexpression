#this script quantile normalizes a directory of files
tic <- as.integer(as.POSIXct(Sys.time()))
library(tidyverse)
library(preprocessCore)
args <- commandArgs(TRUE)
#args[1] is dir to normalize
directory <- args[1]
#create output directory
dirname <- paste0("./", directory, "_quantile_normalized")
if(!dir.exists(dirname)) {
  dir.create(dirname)
}
#list files in dir to normalize
dir_files <-  list.files(directory, full.names = TRUE)

#normalize all files
for (filename in dir_files){
  #read file
  file_data <- read_delim(filename, 
                          delim = "\t", 
                          col_names = TRUE, 
                          col_types = cols(.default = "d", gene = "c")) %>% 
    column_to_rownames("gene") %>% 
    as.matrix()
  #get file basename to name new file
  file_basename <- basename(filename)
  #quantile normalize - normalize.quantiles does not preserve rownames or
  #colnames so have to work around
  normalized_file_data <- normalize.quantiles(file_data)
  #set rownames and colnames 
  colnames(normalized_file_data) <- colnames(file_data)
  rownames(normalized_file_data) <- rownames(file_data)
  #write normalized data to file
  normalized_file_data %>% 
    as.data.frame() %>% 
    rownames_to_column("gene") %>% 
    write_delim(paste0(dirname, "/", 
                       gsub("\\.pcl", "_quantile_norm.pcl", file_basename)),
                delim = "\t",
                col_names = T)
}
#script time
toc <- as.integer(as.POSIXct(Sys.time()))
print(paste("The time it took to run this script in minutes was", (toc-tic)/60, sep = " "))

