# Title: Get Matched Phenotype Data
# Author: Olivia Williamson, Natalie Goulett, Gabriel Odom
# Edited: 2024-04-29



# Install and load packages
library(BiocManager)
library(GEOquery)
library(Biobase)
library(biomaRt)
library(tidyverse)



###  GSE 186332  ###
GSE186332 <- getGEO("GSE186332")
GSE186332_data <- GSE186332[[1]]


head(pData(GSE186332_data))
GSE186332_pheno_names_to_keep <- c(
  "title", "characteristics_ch1.1", "characteristics_ch1.2",
  "characteristics_ch1.3", "treatment_protocol_ch1.1",
  "treatment_protocol_ch1.3", "treatment_protocol_ch1.4",
  "treatment_protocol_ch1.5", "treatment_protocl_ch1.6",
  "experimental set:ch1", "patient response:ch1", "trial arm:ch1"
)
pheno_GSE186332_df <- pData(GSE186332_data) %>%
  select(title, treatment = `characteristics_ch1.1`) %>%
  mutate(SampleID = str_replace(title, pattern = "-", replacement = "\\.")) %>%
  mutate(SampleID = paste0("X", SampleID)) %>%
  select(-title) %>%
  as_tibble()


###  GSE 154041  ###
GSE154041 <- getGEO("GSE154041")
GSE154041_data <- GSE154041[[1]]


head(pData(GSE154041_data))
GSE154041_pheno_names_to_keep <- c(
  "title", "characteristics_ch1", "characteristics_ch1.1", "description",
  "overall survival (months):ch1", "treatment:ch1"
)
pheno_GSE154041_df <- pData(GSE154041_data) %>%
  select(title, treatment = `treatment:ch1`) %>%
  mutate(SampleID = str_extract(title, pattern = "\\d{3}_\\d{2}")) %>%
  mutate(SampleID = paste0("X", SampleID)) %>%
  select(-title) %>%
  as_tibble()

### Merge & Save ###
pheno_bind_df <- bind_rows(pheno_GSE154041_df, pheno_GSE186332_df)
write.csv(
  x = pheno_bind_df,
  file = "./Data_Clean/phenotype_joined_information_20240429.csv",
  row.names = FALSE
)
