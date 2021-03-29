# RNAseq_coexpression
This GitHub repository contains all the code for reproducing the results from our manuscript "**Robust normalization and transformation techniques for constructing gene coexpression networks from RNA-seq data**", which can be found [here](https://doi.org/10.1101/2020.09.22.308577).

The website containing extra results and figures from our analysis can be found [here](https://krishnanlab.github.io/RNAseq_coexpression/index.html).  

# Code
The scripts needed to create a coexpression network with any of our tested workflows can be found in the `src` directory. Most scripts are R scripts and are described below. Some scripts have accompanying data files in the `data` directory that can be used to reproduce our workflow, but can be replaced with custom files to use our workflow with your own data. Calculating the pearson correlations to build the network, the network transformations, and evaluation requires installation of the c++ library, [Sleipnir](https://functionlab.github.io/sleipnir-docs/index.html). The naive and tissue-aware gold standards we used for evaluation can be found in the "gold_standards" directory.

The `docs` and `website` directories are only related to the extra results website (linked above) and nothing should be pulled from these directories for reproduction or use of our workflow.

## 1. Download data, get counts, CPM, RPKM, and TPM
To download the GTEx dataset, use the script `gtex_download_and_normalize.R`. It will read a file from data but requires no arguments. Otherwise, see below.

__Script__:  
`download_and_normalize.R`

__Required R packages__:  
[tidyverse](https://www.tidyverse.org/), [recount](https://bioconductor.org/packages/release/bioc/html/recount.html)

__Purpose__: 
Download entire projects from the [Recount2](https://jhubiostatistics.shinyapps.io/recount/) database and output five directories:
- `counts`: count data for all datasets where the first column in each file is 'gene' and the remaining columns are individual samples 
- `cpm`: cpm-normalized datasets where the first column in each file is 'gene' and the remaining columns are individual samples
- `tpm`: tpm-normalized datasets where the first column in each file is 'gene' and the remaining columns are individual samples
- `rpkm`: rpkm-normalized datasetes where the first column in each file is 'gene' and the remaining columns are individual samples
- `metadata`: available metadata for all datasets
- `rse`: the Rdata object for each dataset
Each remaining step of the pipeline is designed to take a directory as an argument to complete the step for every dataset in the directory, so these directories can now be put through the rest of the workflow as desired.

__Arguments__:  
The single argument is the path to a file that is used to select data for download and normalization. We created this file using the Recount2 metadata that can be accessed with the all_metadata function in the `recount` R package. The metadata was filtered for our desired characteristics and the following columns were selected to write to file: `project`, `sample`, `run`, and `sharq_beta_tissue`.

__Use from command line__:  
```bash
$ Rscript download_and_normalize.R selected_projects-tissue-sample.tsv
```

## 2. Sample selection
We have modified the `download_and_normalize.R` script to only download the data that met our criteria, so this step can be skipped to repeat our analysis. If other data is being used, this step can be executed to filter out samples that have over half zero-expression for lncRNA, antisense RNA, and protein coding genes.

__Script__:  
`sample_filter.R`

__Required R packages__:  
[tidyverse](https://www.tidyverse.org/)

__Purpose__:  
Removes samples from each dataset that have over half zero-expression for lncRNA, antisense RNA, and protein coding genes.

__Arguments__:  
Path to directory of datasets to be sample filtered.

__Use from command line__:  
```bash
$ Rscript sample_filter.R path/directory_to_be_sample_filtered
```

## 3. Gene filtering

__Script__:  
`gene_filtering_by_cpm.R`

__Required R packages__:  
[tidyverse](https://www.tidyverse.org/)

__Purpose__:  
Remove genes which have universally low expression.

__Arguments__:  
The first argument is the path to the directory of datasets to be gene filtered. The second argument is the path to the file that contains the genes to keep. To repeat our analysis with the SRA genes, use `SRA_genes_to_keep_by_cpm_filter.txt` in the `data` directory. To repeat our analysis with the GTEx genes, use `GTEx_genes_to_keep_by_cpm_filter.txt` in the `data` directory. 

__Use from command line__:  
```bash
$ Rscript gene_filtering_by_cpm.R path/directory_to_be_gene_filtered path/genes_to_keep_file.txt
```

## 4. Within-sample normalization
If no within-sample normalization is desired, continue the pipeline with the `counts` directory that has been through any necessary sample selection and gene filtering. Choosing CPM, RPKM, or TPM requires the use of the `cpm`, `rpkm`, or `tpm` directories that have been sample and/or gene filtered as necessary.

## 5. Between-sample normalization
The options are `none`, `TMM`, `upper quartile`, or `quantile` normalization. If no between-sample normalization is desired, skip this step. 

### TMM normalization
TMM normalization requires count data and should not be paired with CPM, RPKM, or TPM.

__Script__:  
`TMM_normalize.R`

__Required R packages__:  
[tidyverse](https://www.tidyverse.org/), [edgeR](https://bioconductor.org/packages/release/bioc/html/edgeR.html)

__Purpose__:  
TMM normalize each dataset in a directory.

__Arguments__:  
Path to directory of datasets to be TMM normalized.

__Use from command line__:  
```bash
$ Rscript TMM_normalize.R path/directory_to_be_normalized
```

### Upper quartile normalization
Upper quartile normalization requires count data and should not be paired with CPM, RPKM, or TPM.

__Script__:  
`upper_quartile_normalize.R`

__Required R packages__:  
[tidyverse](https://www.tidyverse.org/), [edgeR](https://bioconductor.org/packages/release/bioc/html/edgeR.html)

__Purpose__:  
Upper quartile normalize each dataset in a directory.

__Arguments__:  
Path to directory of datasets to be upper quartile normalized.

__Use from command line__:  
```bash
$ Rscript upper_quartile_normalize.R path/directory_to_be_normalized
```

### Quantile normalization
Quantile normalization can be paired with counts, TPM, RPKM, or TPM.

__Script__:  
`quantile_normalize.R`

__Required R packages__:  
[tidyverse](https://www.tidyverse.org/), [preprocessCore](https://bioconductor.org/packages/release/bioc/html/preprocessCore.html)

__Purpose__:  
Quantile normalize each dataset in a directory.

__Arguments__:  
Path to directory of datasets to be quantile normalized.

__Use from command line__:  
```bash
$ Rscript quantile_normalize.R path/directory_to_be_normalized
```

## 6. Gene Type Filtering
This step is not necessary if all gene types are of interest, but the script below will retain only genes types used in our analysis (lncRNA, antisense RNA, and protein coding genes)

__Script__:  
`gene_type_filter_lncrna_antirna_proteincoding.R`

__Required R packages__:  
[tidyverse](https://www.tidyverse.org/)

__Purpose__:  
Remove all gene types that are not lncRNA, antisense RNA, or protein coding.

__Arguments__:  
Path to directory of datasets to be gene type filtered.

__Use from command line__:  
```bash
$ Rscript gene_type_filter_lncrna_antirna_proteincoding.R path/directory_to_be_gene-type_filtered
```

## 7. Hyperbolic arcsine transformation

__Script__:  
`asinh_transform.R`

__Required R packages__:  
[tidyverse](https://www.tidyverse.org/)

__Purpose__:  
Transform data with hyperbolic arcsine function. 

__Arguments__:  
Path to directory of datasets to be transformed.

__Use from command line__:  
```bash
$ Rscript asinh_transform.R path/directory_to_be_transformed
```

## 8. Correlation calculation
This step requires the use of the [Sleipnir](https://functionlab.github.io/sleipnir-docs/index.html) c++ library. Once Sleipnir has been loaded, we use the following command on each individual dataset to output the network edgelist file:

```bash
$ Distancer -i input_dataset.pcl -d pearson -o output_filename.dat -c -s 0 -z
```

This command was incorporated into a simple bash script to iterate over a directory. 

## 9. Network transformation
The options are `none`, `CLR`, or `wTO`. If no network transformation is desired, skip this step.

### CLR
CLR is another option in the [Sleipnir](https://functionlab.github.io/sleipnir-docs/index.html) c++ library. Once Sleipnir has been loaded, we use the following command on each individual dataset to output the CLR transformed network edgelist:

```bash
$ Dat2Dab -i input_dataset.dat -Y -o output_file.dat
```

This command was incorporated into a simple bash script to iterate over a directory. 

### wTO
Our use of the `wTO` R package required that our correlation edgelist be converted to an adjacency matrix before using wTO. For speed, we used a python script to convert the edgelist to an adjacency matrix, then used an R script to use the wTO package and output the transformed edgelist.

#### A - edgelist to adjacency matrix

__Script__:  
`Edgelist_to_matrix.py` (python 3)

__Purpose__:  
Take an edgelist and turn to an adjacency matrix and nodelist for the adjacency matrix.

__Arguments__: 
- dir1 = the directory that the intial edgelist is in  
- edgelist = file in dir1 that should be converted to an adjacency matrix
- dir2 = the directory that the output adjacency matrices and nodelists should be saved to       
- binary = flag that can be used to output a binary rather than weighted adjacency matrix
- diags = what diagonals should be in the adjacency matrix (default = 0)

__Use from command line__:  
```bash
$ python Edgelist_to_matrix.py -dir1 path/input_directory -dir2 path/output_directory -edgelist edgelist_filename.tsv
```

#### B - adjacency matrix to wTO edgelist

__Script__:  
`edgelist_to_wTO_edgelist.R`

__Required R packages__:  
[tidyverse](https://www.tidyverse.org/), [wTO](https://cran.r-project.org/web/packages/wTO/)

__Purpose__:  
wTO tranform the network edges.

__Arguments__:  
The first argument is the path to the directory that output files should be placed in. The second argument is the path to the adjacency matrix output by `Edgelist_to_matrix.py`. The third argument is the path to the nodelist output by `Edgelist_to_matrix.py`.

__Use from command line__:  
```bash
$ Rscript edgelist_to_wTO_edgelist.R path/output_directory adjacency_matrix_file.txt nodelist_file.nodelist
```

## Evaluation
If desired, networks can be evaluated on the functional gold standards described in the manuscript. This requires the use of the [Sleipnir](https://functionlab.github.io/sleipnir-docs/index.html) c++ library. Once Sleipnir has been loaded, we use the following command on each individual dataset:

```bash
$ DChecker -w gold_standard.dab -b 20000 -i input_network.dat -o evaluation_output.txt
```

This command was incorporated into a simple bash script to iterate over a directory. The gold standard can be any gold standard file found in the `gold_standards` directory. The evaluation output will include the number of true positives/false positives/true negatives/false negatives at 20,000 cutoffs so that auPRC, auROC, or other desired metrics can be calculated from the file. In our analysis, the R package [pracma](https://cran.r-project.org/web/packages/pracma/index.html) was used to calculate the area under the precision-recall curve and the area under the ROC curve.

## Data Transformations
VST and rlog are the data transformations we tested to compare to the hyperbolic arcsine transformation. Both transformations can only be used with count data. We perform sample selection and gene filtering by cpm before doing the transformation, then gene type filtering followed by correlation calculation can be performed to build the network. Network transformation can still be done if desired. 

### VST

__Script__:  
`vst_transform.R`

__Required R packages__:  
[tidyverse](https://www.tidyverse.org/), [DESeq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html)

__Purpose__:  
VST transform each dataset in a directory.

__Arguments__:  
Path to directory of datasets to be transformed.

__Use from command line__:  
```bash
$ Rscript vst_transform.R path/directory_to_be_transformed
```

### rlog

__Script__:  
`rlog_transform.R`

__Required R packages__:  
[tidyverse](https://www.tidyverse.org/), [DESeq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html)

__Purpose__:  
rlog transform each dataset in a directory.

__Arguments__:  
Path to directory of datasets to be transformed.

__Use from command line__:  
```bash
$ Rscript rlog_transform.R path/directory_to_be_transformed
```
