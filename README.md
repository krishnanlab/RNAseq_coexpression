# RNAseq_coexpression
This is the GitHub repository for reproducing networks and results from our manuscript "Robust normalization and transformation techniques for constructing gene coexpression networks from RNA-seq data" which can be found here. The website containing extra results and figures from our analysis can be found here. 

# Code
The scripts needed to create a coexpression network with any of our tested workflows can be found in the "src" directory. Most scripts are R scripts and are described below. Some scripts have accompanying data files in the "data" directory that can be used to reproduce our workflow, but can be replaced with custom files to use our workflow with your own data. Calculating the pearson correlations to build the network, the network transformations, and evaluation requires installation of the c++ library, Sleipnir. The naive and tissue-aware gold standards we used for evaluation can be found in the "gold_standards" directory.

The "docs" and "website" directories are only related to the extra results website (linked above) and nothing should be pulled from these directories for reproduction or use of our workflow.

## 1. Download data, get counts, CPM, RPKM, and TPM
__Script__:  download_and_normalize.R

__Required R packages__: tidyverse (CRAN), recount (Bioconductor)

__Purpose__: 
Download entire projects from the Recount2 data base and output five directories:
- counts: count data for all datasets where the first column in each file is 'gene' and the remaining columns are individual samples 
- cpm: cpm-normalized datasets where the first column in each file is 'gene' and the remaining columns are individual samples
- tpm: tpm-normalized datasets where the first column in each file is 'gene' and the remaining columns are individual samples
- rpkm: rpkm-normalized datasetes where the first column in each file is 'gene' and the remaining columns are individual samples
- metadata: available metadata for all datasets
- rse: the Rdata object for each dataset
Each remaining step of the pipeline is designed to take a directory as an argument to complete the step for every dataset in the directory, so these directories can now be put through the rest of the workflow as desired.

__Arguments__:
The single argument is the path to a file that is used to select data for download and normalization. We created this file using the Recount2 metadata that can be accessed with the all_metadata function in the recount R package. The metadata was filtered for our desired characteristics and the following columns were selected to write to file: project, sample, run, sharq_beta_tissue.

__Use from command line__: 
Rscript download_and_normalize.R selected_projects-tissue-sample.tsv

## 2. Sample selection
We have modified the download_and_normalize.R script to only download the data that met our criteria, so this step can be skipped to repeat our analysis. If other data is being used, this step can be executed to filter out samples that have over half zero-expression for lncRNA, antisense RNA, and protein coding genes.

__Script__: sample_filter.R

__Required R packages__: tidyverse (CRAN)

__Purpose__: Removes samples from each dataset that have over half zero-expression for lncRNA, antisense RNA, and protein coding genes.

__Arguments__: Path to directory of datasets to be sample filtered.

__Use from command line__: Rscript sample_filter.R path/directory_to_be_sample_filtered

## 3. Gene filtering

__Script__:  gene_filtering_by_cpm.R

__Required R packages__: tidyverse (CRAN)

__Purpose__: Remove genes which have universally low expression

__Arguments__: The first argument is the path to the directory of datasets to be gene filtered. The decond argument is the path to the file that contains the genes to keep. To repeat our analysis with the SRA genes, use "SRA_genes_to_keep_by_cpm_filter.txt" in the data directory. To repeat our analysis with the GTEx genes, use "GTEx_genes_to_keep_by_cpm_filter.txt" in the data directory. 

__Use from command line__: Rscript gene_filtering_by_cpm.R path/directory_to_be_gene_filtered path/genes_to_keep_file.txt

## 4. Within-sample normalization
If no within-sample normalization is desired, continue the pipeline with the counts directory that has been through any necessary sample selection and gene filtering. Choosing CPM, RPKM, or TPM requires the use of the CPM, RPKM, or TPM directories that have been sample and/or gene filtered as necessary.

## 5. Between-sample normalization
The options are none, TMM, upper quartile, or quantile normalization. If no between-sample normalization is desired, skip this step. 

### TMM normalization
TMM normalization requires count data and should not be paired with CPM, RPKM, or TPM.
__Script__:  

__Required R packages__: tidyverse (CRAN)

__Purpose__: 

__Arguments__:

__Use from command line__: 

### Upper quartile normalization
Upper quartile normalization requires count data and should not be paired with CPM, RPKM, or TPM.
__Script__:  

__Required R packages__: tidyverse (CRAN)

__Purpose__: 

__Arguments__:

__Use from command line__: 

### Quantile normalization
Quantile normalization can be paired with counts, TPM, RPKM, or TPM.
__Script__:  

__Required R packages__: tidyverse (CRAN)

__Purpose__: 

__Arguments__:

__Use from command line__: 

## 6. Gene Type Filtering
This step is not necessary if all gene types are of interest, but the script below will retain only genes types used in our analysis (lncRNA, antisense RNA, and protein coding genes)

__Script__:  

__Required R packages__: tidyverse (CRAN)

__Purpose__: 

__Arguments__:

__Use from command line__: 

## 7. Hyperbolic arcsine transformation
__Script__:  

__Required R packages__: tidyverse (CRAN)

__Purpose__: 

__Arguments__:

__Use from command line__: 

## 8. Correlation calculation
This step requires the use of the [Sleipnir](https://functionlab.github.io/sleipnir-docs/index.html) c++ library. Once Sleipnir has been loaded, we use the following command on each individual dataset to output the network edgelist file:

$ Distancer -i input_dataset.pcl -d pearson -o output_filename.dat -c -s 0 -z

This command was incorporated into a simple bash script to iterate over a directory. 

## 9. Network transformation
The options are none, CLR, or wTO. If no network transformation is desired, skip this step.

### CLR
CLR is another option in the [Sleipnir](https://functionlab.github.io/sleipnir-docs/index.html) c++ library. Once Sleipnir has been loaded, we use the following command on each individual dataset to output the CLR transformed network edgelist:

$ Dat2Dab -i input_dataset.dat -Y -o output_file.dat

This command was incorporated into a simple bash script to iterate over a directory. 

### wTO
Our use of the wTO R package required that our correlation edgelist be converted to an adjacency matrix before using wTO. For speed, we used a python script developed by Anna Yannakopoulos to convert the edgelist to an adjacency matrix, then used an R script to use the wTO package and output the transformed edgelist.
#### A - edgelist to adjacency matrix

__Script__:  

__Required R packages__: tidyverse (CRAN)

__Purpose__: 

__Arguments__:

__Use from command line__: 

#### B - adjacency matrix to wTO edgelist

__Script__:  

__Required R packages__: tidyverse (CRAN)

__Purpose__: 

__Arguments__:

__Use from command line__: 

## Evaluation
If desired, networks can be evaluated on the functional gold standards described in the manuscript. This requires the use of the [Sleipnir](https://functionlab.github.io/sleipnir-docs/index.html) c++ library. Once Sleipnir has been loaded, we use the following command on each individual dataset:

$ DChecker -w gold_standard.dab -b 20000 -i input_network.dat -o evaluation_output.txt

This command was incorporated into a simple bash script to iterate over a directory. The gold standard can be any gold standard file found in the "gold_standards" directory. The evaluation output will include the number of true positives/false positives/true negatives/false negatives at 20000 cutoffs so that auPRC, auROC, or other desired metrics can be calculated from the file.

