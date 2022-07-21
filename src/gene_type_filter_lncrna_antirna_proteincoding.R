#script to filter dirs by only the type of ensembl genes we want
tic <- as.integer(as.POSIXct(Sys.time()))
library(tidyverse)
args <- commandArgs(TRUE)
#args[1] is path to directory needing filtering (also helps name new dir)

#read in ensembl gene ids (1st col) and type of gene (2nd col) 
gene_info <- read_delim("../data/ensembl-id_gene-type_mygeneinfo.tsv",
                        delim = "\t",
                        col_names = T)
#filter for types of genes we want
selected_gene_info <- gene_info %>% 
  filter(type_of_gene %in% c("protein_coding", "antisense_RNA", "lincRNA"))
#pull ensembl ids only
selected_ensembl_ids <- selected_gene_info %>% pull(ensembl_ids)

#dir files to be filtered
dir_files <- list.files(args[1], full.names = TRUE)
#create dir to write files to
dirname <- paste0("./", base::basename(args[1]), "_gene_filtered")
if(!dir.exists(dirname)) {
  dir.create(dirname)
}

#filter all files
for (file_name in dir_files){
  file_data <- read_delim(file_name, 
                          delim = "\t", 
                          col_names = TRUE, 
                          col_types = cols(.default = "d", gene = "c"))
  #get basename
  file_basename <- base::basename(file_name)
  #get rid of decimal points in ensembl ids in files (pull column, split on ".")
  ensembl_id_with_dec <- pull(file_data, var = 1)
  ensembl_id_dec_removed <- vapply(strsplit(ensembl_id_with_dec,".",fixed = TRUE), `[`, 1, FUN.VALUE=character(1))
  #remove gene column from data file to replace with gene column w/o decimal
  file_data_no_gene_col <- file_data %>% select(-1)
  #put ensembl id col w/o decimal at start of df, write to file
  no_dec_ensembl_ids_df <- data.frame(gene=ensembl_id_dec_removed, file_data_no_gene_col)
  no_dec_ensembl_ids_df %>% 
    as_tibble() %>% 
    filter(gene %in% selected_ensembl_ids) %>% 
    as.data.frame() %>% 
    write_delim(paste0(dirname, "/", gsub("\\.pcl$", "_gene-type_filtered.pcl", file_basename)),
                delim = "\t",
                col_names = TRUE)
}
#script time
toc <- as.integer(as.POSIXct(Sys.time()))
print(paste("The time it took to run this script in minutes was", (toc-tic)/60, sep = " "))

