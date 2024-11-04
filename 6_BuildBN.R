# Title: Bayesian Network Analysis with bnlearn
# Author: Olivia Williamson, Natalie Goulett, Gabriel Odom
# Created: 2024-04-26
# Edited: 2024-04-28



# Install and load packages
# install.packages("BiocManager")
# BiocManager::install("DESeq2")
# install.packages("openxlsx")
# install.packages("bnlearn")
# BiocManager::install("Rgraphviz")
library(bnlearn)
library(Rgraphviz)
library(tidyverse)

# Grabbing file from previous script - "Merge_Remove_Batch_Effects"
final_data <- readRDS(file = "Data_Clean/GSE_discretized_20240426.rds")
# final_data_fct <- as.factor(final_data)



##########  join to treatment data #############################################

pheno_df <- read.csv("Data_Clean/phenotype_joined_information_20240429.csv")

# Assuming pheno_df already has the 'treatment' column loaded
# Convert the 'treatment' column to a set of binary columns for each treatment group
# pheno_df_treatments <- pheno_df %>%
#   mutate(across(treatment, ~ factor(
#     .x,
#     levels = c(
#       "Lomustine",
#       "Regorafenib",
#       "trial arm: A",
#       "trial arm: B",
#       "trial arm: C",
#       "trial arm: D"
#     )
#   ))) %>%
#   mutate(across(treatment, as.integer)) %>%
#   pivot_wider(names_from = treatment, values_from = treatment, values_fill = 0)

pheno_df_treatments <- pheno_df %>%
  mutate(
    Lomustine = as.integer(treatment == "Lomustine"),
    Regorafenib = as.integer(treatment == "Regorafenib"),
    `Selinexor Pre-Op` = as.integer(treatment == "trial arm: A"),
    `Selinexor 50mg` = as.integer(treatment == "trial arm: B"),
    `Selinexor 60mg` = as.integer(treatment == "trial arm: C"),
    `Selinexor 80mg` = as.integer(treatment == "trial arm: D")
  ) %>%
  select(-treatment)  # Remove the original treatment column if it's no longer needed

# Assuming that the new treatment columns have been correctly added, you may not need to rename them as the above code assigns the names directly.


# Renaming columns
new_treatment_names <- c("Lomustine", "Regorafenib", "Selinexor Pre-Op", "Selinexor 50mg", "Selinexor 60mg", "Selinexor 80mg")
names(pheno_df_treatments)[2:7] <- new_treatment_names


# test with sample
# final_data2 <- as.data.frame(final_data)
# finaldata_sample <- sample(final_data2, 10)

#set.seed(8888)
# create random genes for testing code
# random_genes_idx <- sample(x = 1:nrow(final_data), size = 10)
gbm_pathway_lit <- c("POT1", "HERC2", "BRIP1", "POLE", "EGFR", "SOX9", "PGK1",
                     "CA9", "VEGFA", "SPP1", "HIF1A", "HP1BP3", "ZC3H7B", "BRAF",
                     "PTEN", "MGMT", "PMS2", "IDH1", "NF1", "TP53", "PDGFRA",
                     "AEBP1", "ANXA2R", "TMEM60", "PRRG3", "RPS4X", "CTSK",
                     "SLC2A1", "KLHL12", "CDKN1A", "CA12", "WDR1", "CD53",
                     "CBR4", "NIFK-AS1", "RAB30-DT", "FOXO1", "NGFR")

final_data_sample <- as.data.frame(
  # bnlearn requires a data frame in sample x gene format
  t(
    final_data[gbm_pathway_lit,]
  )
)

# Join the new treatment columns to the data
final_joined_df <- final_data_sample %>%
  tibble::rownames_to_column("SampleID") %>%
  left_join(pheno_df_treatments, ., by = "SampleID") %>%
  tibble::column_to_rownames(var = "SampleID")

# final_joined_df <- final_data_sample %>%
#   tibble::rownames_to_column("SampleID") %>%
#   left_join(pheno_df, ., by = "SampleID") %>%
#   tibble::column_to_rownames(var = "SampleID")

# keep track of sample ids
sample.id <- rownames(final_joined_df)
final_data_fct <- data.frame(
  lapply(final_joined_df,
         FUN = function(x) factor(x))
)
str(final_data_fct)

# re-add sample ids to row names
rownames(final_data_fct) <- sample.id
rm(sample.id)

#diff exp analysis/around 200 goi/ from GBM literature
a <- Sys.time()
gene_hc <- hc(final_data_fct, restart = 5000)
b <- Sys.time()

graphviz.plot(gene_hc, layout = "dot")
gene_MLE <- bn.fit(gene_hc, data = final_data_fct, method = "mle")
gene_MLE

