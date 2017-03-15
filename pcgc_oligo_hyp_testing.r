# hypothesis testing for significance of oligogenic hypothesis in PCGC trios
# Kathleen Keough 2.16.17

# load necessary packages
setwd('/Users/kathleen/Box Sync/UCSF/ConklinPollard/projects/PCGC_oligo_proj')
library(dplyr)
library(readr)

# this is your file with 1/0 hypotheses for genes of interest (or all genes!)

trio_dat <- read_tsv("cardio277_table_freq30_excNAs.txt")
print(colnames(trio_dat[,apply(trio_dat==0,2,all)]))
print(dim(trio_dat))
# filter genes where all zeros indicating no trios met oligo hypothesis

trio_dat <- trio_dat[ , !apply(trio_dat==0,2,all)]
print(dim(trio_dat))
# construct gene list from column names (column 3 onwards)

list_of_interest = colnames(trio_dat) 
gene_list <- as.list(tail(list_of_interest,length(list_of_interest)-2))

# perform Fisher's exact test for each gene

results_df <- data.frame(matrix(ncol = 2, nrow = length(gene_list)))
colnames(results_df) <- c("Genes","Fisher's Exact p-value")
results_df$Genes <- gene_list
rownames(results_df) <- results_df$Genes

for (gene in gene_list){
  fishers_pval <- fisher.test(trio_dat$CaseControl,
                              y=trio_dat[[paste(gene)]],
                              alternative = "greater")$p.value
  results_df[gene,2] <- fishers_pval
}

results_df$Genes <- unlist(results_df$Genes)

write.table(results_df, file='277cardio_common+rare_results.tsv', sep='\t', row.names = FALSE)
