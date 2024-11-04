# Title: Ensure both datasets have matchable gene symbols
# Author: Olivia Williamson, Natalie Goulett, Gabriel Odom
# Created: 2024-04-10
# Edited: 2024-04-28



# install and load the following packages

# install.packages("hgnc")
# install.packages("tidyverse")
# BiocManager::install("edgeR")
library("hgnc")
library("tidyverse")
library("edgeR")


# import gene name crosswalk table containing the four gene name formats
hgnc_tbl <- import_hgnc_dataset() %>%
  select("entrez_id", "ensembl_gene_id", "hgnc_id", "symbol")

# GSE_154041_long
# Reads counts were normalized with the TMM normalization method and
# subsequently for reads depth as count per million (CPM) using
# the cpm function implemented in edgeR. (**what does this mean?)

# read in GSE154041 data from the DownloadData R script
GSE154041_long_df <- read.delim(
  "GSE154041/GSE154041_samples_expression_table_not_filtered.txt.gz"
)
# create GSE153041 dataset with gene symbols as name
GSE154041Final <- GSE154041_long_df %>%
  rename("ensembl_gene_id" = "Ensembl_geneid") %>%
  select(-"Gene_description", -"Gene_type") %>%
  select("Gene_name", everything())
sum(is.na(GSE154041Final$Gene_name))

# **what was this anti_join for? To test for any genes that don't match?
#   Shall we keep it?
# anti_join(GSE154041Final, hgnc_tbl, by = "ensembl_gene_id") %>%
#   head()
write_csv(GSE154041Final, file = "Data_Clean/GSE154041_wsymbols_20240410.csv")


# GSE_186332_long
# read in GSE154041 data from the DownloadData R script
GSE186332_long_df <- read.delim(
  "GSE186332/GSE186332_gbm_rawcounts.tsv.gz"
)
GSE186332Final <- GSE186332_long_df %>%
  rename("entrez_id" = "NCBIGeneID") %>%
  # Induces 17 NAs, problems with ERCC-00002, ..., ERCC-00033; this is fine: we
  #   checked and these NA genes are all technical controls
  mutate(entrez_id = as.integer(entrez_id)) %>%
  filter(!is.na(entrez_id))

# create GSE186332 dataset with gene symbols as name
GSE186332_symbols <-
  hgnc_tbl %>%
  left_join(GSE186332Final, ., by = "entrez_id") %>%
  select("symbol", "ensembl_gene_id", "entrez_id", "hgnc_id", everything()) %>%
  rename("Gene_name" = "symbol") %>%
  select(-entrez_id, -hgnc_id)

write_csv(
  GSE186332_symbols,
  file = "Data_Clean/GSE186332_wsymbols_20240410.csv"
)