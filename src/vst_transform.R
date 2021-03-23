tic <- as.integer(as.POSIXct(Sys.time()))
library(tidyverse)
library(DESeq2)
args <- commandArgs(TRUE)
#args[1] is directory of files to transform

directory <- args[1]
#create output dir
output_dir <- paste0("./", directory, "_VST_transformed")
if(!dir.exists(output_dir)) {
  dir.create(output_dir)
}
#list files in dir to normalize
dir_files <- list.files(directory, full.names = TRUE)

for (filename in dir_files){
  #read in file, convert to matrix for norm factor function
  count_data <- read_delim(filename, 
                           delim = "\t", 
                           col_names = T, 
                           col_types = cols(.default = "d", gene = "c")) %>% 
    column_to_rownames("gene") %>% 
    as.matrix()
  #get file basename
  file_basename <- basename(filename)
  #vst transform
  vst_trans <- vst(count_data)
  #write to file
  vst_trans %>% 
    as.data.frame() %>% 
    rownames_to_column("gene") %>% 
    write_delim(paste0(output_dir, "/", gsub("\\.pcl$", "_VST_trans.pcl", file_basename)),
                delim = "\t",
                col_names = T)
}


#script time
toc <- as.integer(as.POSIXct(Sys.time()))
print(paste("The time it took to run this script in minutes was", (toc-tic)/60, sep = " "))
