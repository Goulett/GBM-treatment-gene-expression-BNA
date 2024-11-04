# Title: Download Data from GEO
# Author: Olivia Williamson & Natalie Goulett
# Created: 2024-04-04
# Edited: 2024-04-28



#We will acquire data from two NCBI datasets for our Bayesian network analysis -
#   GSE186332 and GSE154041.
# GSE154041's genes are named using ensemble IDs whereas GSE186332 contains NCBI
#   Gene names. We will convert them both to Hugo gene symbols in the following 
#   R script (step 2).
# GSE154041 read counts were normalized with the TMM normalization method and
#   subsequently filtered(?)** for reads depth as count per million (CPM) using 
#   the cpm function implemented in edgeR. GSE186332 contains raw count data, so
#   we will normalize it in the normalization R script (step 3).
# More info on the data are available at 
#   https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE154041 and
#   https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi.

# Install and load the following packages

# install.packages("BiocManager")
# BiocManager::install("GEOQuery")
# BiocManager::install("Biobase")
# BiocManager::install("biomaRt")
# install.packages("tidyverse")
# install.packages("hgnc")
library(BiocManager)
library(GEOquery)
library(Biobase)
library(biomaRt)
library(tidyverse)
library("hgnc")


########################## READING IN RAW COUNTS ###############################

####### GSE186332
# Gene features (obs. count) - 25,132 obs
# Samples - 57
# treatment column - "GSE186332_pheno_names_to_keep$treatment_protocol_ch1.1"

# Treatment - Selinexor doses (twice weekly) preoperatively (Arm A; n = 8
#   patients). Patients not undergoing surgery received 50 mg/m2 (Arm B, n = 24),
#   or 60 mg (Arm C, n = 14) twice weekly, or 80 mg once weekly (Arm D; n = 30).
# Phenotype info available - tissue type; trial arm; experimental set
#   (Discovery, validation); patient response (responder, non-responder);
#   treatment protocol (keep 1.1-1.6) including drug type, dose, frequency;
#   Primary endpoint was 6-month progression-free survival rate (PFS6).
# Gene name format: NCBI ID, which is Entrez ID. Convert to Symbol format using:
# https://support.bioconductor.org/p/125380/
# Reference Genome: hg38
GSE186332 <- getGEO("GSE186332")
GSE186332_data <- GSE186332[[1]]
# gets formal class ExpressionSet`
GSE186332_data
names(pData(GSE186332_data))
head(pData(GSE186332_data))
GSE186332_pheno_names_to_keep <- c(
  "title", "characteristics_ch1.1", "characteristics_ch1.2",
  "characteristics_ch1.3", "treatment_protocol_ch1.1",
  "treatment_protocol_ch1.3", "treatment_protocol_ch1.4",
  "treatment_protocol_ch1.5", "treatment_protocl_ch1.6",
  "experimental set:ch1", "patient response:ch1", "trial arm:ch1"
)
GSE186332_raw <- getGEOSuppFiles("GSE186332")
# data was downloaded to top level directory. It was manually moved to Data/.
GSE186332_raw
GSE186332_long_df <- read.delim(
  "GSE186332/GSE186332_gbm_rawcounts.tsv.gz"
)
# use version GPL24676, which is the datafile for GSE186332.
head(GSE186332_long_df)

# **create a vector of ensemble IDs from homo sapiens set and then ___?
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
mapping <- getBM(attributes = c())

# Genes from the GSE186332 dataset are named using NCBI gene IDs, which are
#   identical to Entrez IDs. Here we use the `AnnotationHub` and `HGNC` packages
#   to add Hugo gene symbols to the dataset.
import_hgnc_dataset()
head(tbl)
GSE186332Final <- GSE186332_long_df %>%
  rename("entrez_id" = "NCBIGeneID")
hgnc_tbl <- import_hgnc_dataset() %>%
  select("entrez_id", "ensembl_gene_id", "hgnc_id", "symbol")
hgnc_tbl %>% as.data.frame()
didthiswork <- merge(hgnc_tbl, GSE186332Final, by = "entrez_id")


###### GSE154041
# 74 Samples - 10902 observations (filtered)
# treatment column - "GSE154041_pheno_names_to_keep$treatment:ch1"

# Treatments - Regorafenib or Lomustine
# Phenotype info available - Patient ID, treatment type, and survival (in 
#   months).
# Reference Genome - HG 38
# Data are have been normalized using TMM method.
GSE154041 <- getGEO("GSE154041")
class(GSE154041)
length(GSE154041)
names(GSE154041)
GSE154041_data <- GSE154041[[1]]

names(pData(GSE154041_data))
head(pData(GSE154041_data))
names(pData(GSE154041))
GSE154041_pheno_names_to_keep <- c(
  "title", "characteristics_ch1", "characteristics_ch1.1", "description",
  "overall survival (months):ch1"
)
GSE154041_raw <- getGEOSuppFiles("GSE154041")
GSE154041_raw
GSE154041_long_df <- read.delim(
  "GSE154041/GSE154041_samples_expression_table_filtered.txt.gz"
)
head(GSE154041_long_df)

# GSE154041 genes are named using ensemble IDs whereas GSE186332 contains NCBI
#   Gene names. 
# rename the ensemble ID column to match format of hgnc table(?)**
GSE154041Final <- GSE154041_long_df %>% 
  rename("ensembl_gene_id" = "Ensembl_geneid")
# Add hgnc table names to GSE154041
didthiswork2 <- merge(
  hgnc_tbl,
  GSE154041Final,
  by = "ensembl_gene_id",
  .before = "TRUE"
)
