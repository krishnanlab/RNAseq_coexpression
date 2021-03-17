The dataframes necessary for each Rmd page are in this directory.
They were all created by using the data preprecessing procedures in
/Users/kayla/rnaseq_coexp/results/2020-02-04_RNAseq_coexpression_paper_figure_1.Rmd
before line 478, where it says ######END DATA PREPROCESSING###############
Dataframes saved:
* metadata (all metadata of projects/samples)
* method_codes 
* gtex_naive (GTEx naive analysis - log2_auprc_over_prior, auroc)
* gtex_knowledge (GTEx tissue-specfic analysis - log2_auprc_over_prior, auroc)
* sra_naive (SRA naive analysis - log2_auprc_over_prior, auroc)
* sra_knowledge (SRA tissue-specfic analysis - log2_auprc_over_prior, auroc)
* grnaive (GTEx resampling naive analysis - log2_auprc_over_prior, auroc)
* grknowledge (GTEx resampling tissue-specfic analysis - log2_auprc_over_prior, auroc)
* prnaive (GTEx and SRA prec-recall naive results)
* prknowledge (GTEx and SRA prec-recall tissue-specific results)
* prgr (GTEx and SRA prec-recall naive results)
* prgrk (GTEx and SRA prec-recall tissue-specific results)

prnaive and prknowledge were written as is but also separated into GTEx (2020-03-19_prec_recall_gtex_[naive|knowledge]_results.txt) 
and SRA (2020-03-19_prec_recall_sra_[naive|knowledge]_results.txt) results (filter by project, write to file)
