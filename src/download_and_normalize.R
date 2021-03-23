tic <- as.integer(as.POSIXct(Sys.time()))

library("tidyverse")
library("recount")
#argument should be file made from recount2 metadata (see README)
args <- commandArgs(TRUE)

# make some dirs if they don't exist
#####################################
dirname <- "./rse"
if(!dir.exists(dirname)) {
  dir.create(dirname)
}

dirname <- "./counts"
if(!dir.exists(dirname)) {
  dir.create(dirname)
}

dirname <- "./metadata"
if(!dir.exists(dirname)) {
  dir.create(dirname)
}

dirname <- "./rpkm"
if(!dir.exists(dirname)) {
  dir.create(dirname)
}

dirname <- "./cpm"
if(!dir.exists(dirname)) {
  dir.create(dirname)
}

dirname <- "./tpm"
if(!dir.exists(dirname)) {
  dir.create(dirname)
}

# Functions for calulating counts, cpm, rpkm, and tpm from recount2 rse
########################################################################

#Get counts data_frame from recount2 rse
get_counts <- function(rse){
  libtype_factor <- rep(1, length(colData(rse)$paired_end))
  libtype_factor[colData(rse)$paired_end] <- 2
  
  round(sweep(assays(rse)$counts, 2, libtype_factor*colData(rse)$avg_read_length, "/")) %>%
    as.data.frame()
}


#Calculate cpm from recount2 rse
get_cpm <- function(tis_gene_count) {
  sweep(tis_gene_count, 2, c(colSums(tis_gene_count)/(10**6)), "/") %>%
    apply(., 2, round, digits = 5) %>%
    as.data.frame()
}


#Calculate tpm from recount2 rse 
get_tpm <- function(tis_gene_count){
  rpk <- sweep(tis_gene_count, 1, c((bplength/1000)), "/")
  sweep(rpk, 2, c((colSums(rpk)/(10**6))), "/") %>%
    apply(., 2, round, digits = 5) %>%
    as.data.frame()
}


#Calculate rpkm from recount2 rse
get_rpkm <- function(tis_gene_count){
  sweep(get_cpm(tis_gene_count), 1, c((bplength/1000)), "/") %>%
    apply(., 2, round, digits = 5) %>%
    as.data.frame()
}

# Get selected data
#######################
selected_projects_sample_tissue <- read_delim(args[1],
                                              delim = "\t",
                                              col_names = T)
#make unique projects = selected_projects vector
selected_projects <- selected_projects_sample_tissue %>%
  pull(project) %>%
  unique()

for (proj in selected_projects){
  print(proj)
  download_study(proj, type = "rse-gene")
  file.rename(paste0(proj,"/rse_gene.Rdata"),
              paste0("./rse/", proj, "_rse_gene.Rdata"))
  
  load(file.path("./rse", paste0(proj, "_rse_gene.Rdata")))
  
  #keep all metadata for rse, file name = "[projID]_metadata.txt"
  rse_gene %>%
    colData() %>%
    as.data.frame() %>%
    write_delim(paste0("./metadata/", paste(proj, "metadata.txt", sep = "_")),
                delim = "\t",
                col_names = T)
  # Extract bp_length from rowranges
  bplength <- rse_gene %>%
    rowRanges() %>%
    as_tibble() %>%
    pull(bp_length)
  
  # Get gene_counts per tissue and normalize the counts matrix
  ############################################################
  #convert rse coverage to actual gene counts with get_counts()
  gene_counts <- rse_gene %>%
    get_counts()
  
  #pull unique tissues from the project 
  tissues <- selected_projects_sample_tissue %>%
    filter(project == proj) %>%
    pull(sharq_beta_tissue) %>%
    unique()
  
  #filter for the samples (runs) we actually want from each project (using project and desired tissue from that project)
  for (tis in tissues) {
    # Get samples corresponding to project and tissue
    runs <- selected_projects_sample_tissue %>%
      filter(project == proj & sharq_beta_tissue == tis) %>%
      pull(run)
    
    
    # Get gene_counts for given project and tissue
    gene_counts %>% 
      select(intersect(runs, colnames(gene_counts))) %>% 
      rownames_to_column(var = "gene") %>%
      write_delim(paste0("./counts/", paste(proj, gsub(" ", "_", tis), "counts.txt", sep = "_")),
                  delim = "\t",
                  col_names = T)
    
    
    #tissue_gene_count = argument for get_cpm, get_tpm, and get_rpkm
    tissue_gene_count <- read_delim(paste0("./counts/", paste(proj, gsub(" ", "_", tis), "counts.txt", sep = "_")),
                                    delim = "\t",
                                    col_names = T,
                                    col_types = cols(.default = "d", gene = "c")) %>%
      column_to_rownames("gene")
    
    
    #get cpm, write to file
    get_cpm(tissue_gene_count) %>%
      rownames_to_column(var = "gene") %>%
      write_delim(paste0("./cpm/", paste(proj, gsub(" ", "_", tis), "cpm.txt", sep = "_")), 
                  delim = "\t", 
                  col_names = T)
    
    
    #get rpkm, write to file
    get_rpkm(tissue_gene_count) %>%
      rownames_to_column(var = "gene") %>%
      write_delim(paste0("./rpkm/", paste(proj, gsub(" ", "_", tis), "rpkm.txt", sep = "_")), 
                  delim = "\t", 
                  col_names = T)
    
    
    #get tpm, write to file
    get_tpm(tissue_gene_count) %>%
      rownames_to_column(var = "gene") %>%
      write_delim(paste0("./tpm/", paste(proj, gsub(" ", "_", tis), "tpm.txt", sep = "_")), 
                  delim = "\t", 
                  col_names = T)
  }
}


toc <- as.integer(as.POSIXct(Sys.time()))
print(paste("The time it took to run this script in minutes was", (toc-tic)/60, sep = " "))

