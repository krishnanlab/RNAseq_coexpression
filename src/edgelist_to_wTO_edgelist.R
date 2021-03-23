tic <- as.integer(as.POSIXct(Sys.time()))

library(tidyverse)
library(wTO)

args <- commandArgs(TRUE)
directory <- args[1]
adjfile <- args[2]
nodefile <- args[3]

file_basename <- basename(adjfile)
print("adjfile:")
print(adjfile)
print("nodefile:")
print(nodefile)
#read in file
adj_matrix  <- read_delim(adjfile, delim = " ", col_names = F) %>% as.matrix()
node_list <- read_delim(nodefile, delim = "\t", col_names = F) %>% pull(1)
#read file time
toc <- as.integer(as.POSIXct(Sys.time()))
print(paste("The time it took to read both files in minutes was", (toc-tic)/60, sep = " "))

#give colname and rownames to matrix
colnames(adj_matrix) <- node_list
rownames(adj_matrix) <- node_list

wto_matrix <- wTO(adj_matrix, sign = "sign")
#wTO time
toc <- as.integer(as.POSIXct(Sys.time()))
print(paste("The time it took to finish reading and wTO in minutes was", (toc-tic)/60, sep = " "))

wto_edge_list <- wTO.in.line(wto_matrix)

#write new egde list to file
wto_edge_list %>% 
  write_delim(paste0(directory, "/", gsub("\\.dat_AdjMatrix_Weighted_Diags0.00.txt", "_wTO_edgelist.dat", file_basename)),
              delim = "\t", col_names = F)
#script time
toc <- as.integer(as.POSIXct(Sys.time()))
print(paste("The time it took to run this script in minutes was", (toc-tic)/60, sep = " "))
