# RNAseq_coexpression
This is the GitHub repository for reproducing networks and results from our manuscript "Robust normalization and transformation techniques for constructing gene coexpression networks from RNA-seq data" which can be found here. The website containing extra results and figures from our analysis can be found here. 

# Code
All of the scripts needed to create a coexpression network with any of our tested workflows can be found in the "source_code" directory. Most scripts are R scripts and are described below. Some scripts have accompanying data files in the "data" directory that can be used to reproduce our workflow, or can be replaced with custom files to use our workflow with your own data. Calculating the pearson correlations to build the network, the network transformations, and evaluation requires installation of the c++ library, Sleipnir. The naive and tissue-aware gold standards we used for evaluation can be found in the "gold_standards" directory.

## 1. Download data, get counts, CPM, RPKM, and TPM
__Script__:  download_and_normalize.R
__Required R packages__: tidyverse (CRAN), recount (Bioconductor)
__Purpose__: 
Download entire projects from the Recount2 data base and output five directories:
- counts: count data for all datasets where the first column is 'gene' and the remaining columns are individual samples
- cpm: cpm-normalized datasets where the first column is 'gene' and the remaining columns are individual samples
- tpm: tpm-normalized datasets where the first column is 'gene' and the remaining columns are individual samples
- rpkm: rpkm-normalized datasetes where the first column is 'gene' and the remaining columns are individual samples
- metadata: available metadata for all datasets
Each remaining step of the pipeline is designed to take a directory as an argument to complete the step for every dataset in the directory, so these directories can now be put through the rest of the workflow as desired.
__Arguments__:
The single argument is a file that is used to select data for download and normalization. It has the following columns: project ID, tissue, sample ID. All metadata from Recount2 can be downloaded and filtered to create this file if using a different set of data than our provided file (see all_metadata() function from recount R package). 
__Use from command line__: 
Rscript download_and_normalize.R selected_projects-tissue-sample.tsv

## 2. Sample selection
__Script__:  
__Required R packages__: tidyverse (CRAN)
__Purpose__: 
__Arguments__:
__Use from command line__: 

## 3. Gene filtering
__Script__:  
__Required R packages__: tidyverse (CRAN)
__Purpose__: 
__Arguments__:
__Use from command line__: 

## 4. Within-sample normalization
If no within-sample normalization is desired, continue the pipeline with the counts directory that has been through any necessary sample selection and gene filtering. Choosing CPM, RPKM, or TPM requires the use of the CPM, RPKM, or TPM directories that have been sample and/or gene filtered as necessary.

## 5. Between-sample normalization
The options are none, TMM, upper quartile, or quantile normalization. If no between-sample normalization is desired, skip this step. 
### A - TMM normalization
__Script__:  
__Required R packages__: tidyverse (CRAN)
__Purpose__: 
__Arguments__:
__Use from command line__: 

### B - upper quartile normalization
__Script__:  
__Required R packages__: tidyverse (CRAN)
__Purpose__: 
__Arguments__:
__Use from command line__: 

### C - quantile normalization
__Script__:  
__Required R packages__: tidyverse (CRAN)
__Purpose__: 
__Arguments__:
__Use from command line__: 

