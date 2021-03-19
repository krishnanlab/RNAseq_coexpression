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

## Correlation calculation
This step requires the use of the [Sleipnir](https://functionlab.github.io/sleipnir-docs/index.html) c++ library. Once Sleipnir has been loaded, we use the following command on each individual dataset to output the network edgelist file:
$ Distancer -i input_dataset.pcl -d pearson -o output_filename.dat -c -s 0 -z
This command was incorporated into a simple bash script to iterate over a directory. 

## Network transformation
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

