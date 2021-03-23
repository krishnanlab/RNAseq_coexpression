tic <- as.integer(as.POSIXct(Sys.time()))
library(tidyverse)
args <- commandArgs(TRUE)
#args[1] is path to directory that must be sample filtered

file_vector <- list.files(path = args[1], full.names = TRUE)
#read in first file 
count_values <- read_delim(file_vector[1],
                           delim = "\t",
                           col_names = T,
                           col_types = cols(.default = "d", gene = "c"))

#read in remaining files, select all but gene column, cbind to large data frame
for (file_name in file_vector[-1]){
  new_count_values <- read_delim(file_name, 
                                 delim = "\t",
                                 col_names = T,
                                 col_types = cols(.default = "d", gene = "c")) %>% 
    select(-1)
  #cbind every file to make large df
  count_values <- cbind(count_values, new_count_values)
}

#get rid of decimal points in ensembl ids in files (pull column, split on ".", select first part)
ensembl_id_with_dec <- pull(count_values, var = 1)
ensembl_id_dec_removed <- vapply(strsplit(ensembl_id_with_dec,".",fixed = TRUE), `[`, 1, FUN.VALUE=character(1))
#add column without decimal point to tibble with all samples ("dec_removed")
count_values <- cbind(dec_removed = ensembl_id_dec_removed, count_values)

#read in ensembl gene ids (1st col) and type of gene (2nd col) 
gene_info <- read_delim("/mnt/home/john3491/projects/rnaseq_coexp/data/ensembl-id_gene-type_mygeneinfo.txt",
                        delim = "\t",
                        col_names = T)
#filter for types of genes we want
selected_gene_info <- gene_info %>% 
  filter(type_of_gene %in% c("protein_coding", "antisense_RNA", "lincRNA"))
#pull ensembl ids only
selected_ensembl_ids <- selected_gene_info %>% pull(ensembl_ids)

#filter sample df to only genes we want to look at, get rid of gene columns
filtering_df <- count_values %>% 
  filter(dec_removed %in% selected_ensembl_ids)
#find what number half of genes of interest contained in projects is
cutoff <- length(unique(filtering_df$dec_removed))/2
#get rid of gene columns (first 2) in filtering df so it's numbers only
filtering_df <- filtering_df %>% select(-1) %>% select(-1)

#find number of zero-exp genes  of interest in each sample
zero_vector <- colSums(filtering_df == 0)
#find samples with over half zero-exp genes of interest
samples_to_remove <- zero_vector[zero_vector > cutoff]

#create output dir if it doesn't exist
directory_basename <- basename(args[1])
dirname <- paste0("./", directory_basename, "_sample_filtered")
if(!dir.exists(dirname)) {
  dir.create(dirname)
}

for (filename in file_vector){
  file_data <- read_delim(filename, 
                          delim = "\t", 
                          col_names = TRUE, 
                          col_types = cols(.default = "d", gene = "c"))
  #get basename for writing new file
  file_basename <- base::basename(filename)
  #select only desired samples (runs)
  selected_data <- select(file_data, -one_of(samples_to_remove))
  #write file again without 
  selected_data %>% 
    as_data_frame() %>% 
    write_delim(paste0(dirname, "/", gsub("\\.txt$", "_sample_filtered.pcl", file_basename)),
                delim = "\t",
                col_names = TRUE)
}

#script time
toc <- as.integer(as.POSIXct(Sys.time()))
print(paste("The time it took to run this script in minutes was", (toc-tic)/60, sep = " "))
