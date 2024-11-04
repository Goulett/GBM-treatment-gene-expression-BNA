# Bayesian Network Analysis of Treatment-Gene Expression Relationships in Glioblastoma Multiforme

## Authors

- [Olivia Williamson](https://github.com/oliviawilliamson)  
- [Natalie Goulett](https://github.com/Goulett)  
- [Gabriel Odom](https://github.com/gabrielodom)  

## Overview

This repository contains R scripts and associated data for a research project analyzing gene expression and treatment effects in patients with Glioblastoma Multiforme (GBM). Using two datasets from the Gene Expression Omnibus (GEO) database, we constructed a Bayesian network to explore relationships between three chemotherapy treatments (Regorafenib, Lomustine, and Selinexor) and gene expression within these patients.

The project pipeline, as outlined in *Comparing Three Chemotherapy Treatments' Effects on Gene Expression in Patients with Glioblastoma Multiforme using the Hill Climbing Algorithm*, by Williamson, O. and Goulett, N. (2024), involves data acquisition, preprocessing, normalization, merging, batch effect removal, and Bayesian network analysis. The Bayesian network was built using the `bnlearn` package with the Hill Climbing algorithm.

## Scripts

Each R script corresponds to a specific step in the analysis pipeline. The scripts are designed to be run sequentially.

### 1. `1_DownloadData.R`
- **Purpose**: Download raw count data from GEO for two datasets: GSE186332 (raw counts, normalized later) and GSE154041 (pre-normalized with TMM).
- **Authors**: Olivia Williamson & Natalie Goulett
- **Details**: Raw counts and metadata are downloaded and stored locally. Gene names are converted to Hugo gene symbols for consistency across datasets.

### 2. `2_AddGeneSymbols.R`
- **Purpose**: Convert gene identifiers to consistent Hugo gene symbols across both datasets.
- **Authors**: Olivia Williamson, Natalie Goulett, Gabriel Odom
- **Details**: Utilizes the `hgnc` package to map gene symbols, preparing datasets for subsequent analysis.

### 3. `3_NormalizationDiscretization.R`
- **Purpose**: Normalize gene expression data and discretize it into categories of over-, non-, and under-expressed genes.
- **Authors**: Olivia Williamson, Natalie Goulett, Gabriel Odom
- **Details**: The TMM method is applied to raw counts, then counts per million (CPM) are calculated and scaled. Data is discretized based on thresholds.

### 4. `4_MergeRemoveBatchEffects.R`
- **Purpose**: Merge datasets and remove batch effects to ensure compatibility for network analysis.
- **Authors**: Olivia Williamson, Natalie Goulett, Gabriel Odom
- **Details**: Uses `limma` to remove batch effects between GSE154041 and GSE186332, ensuring harmonized data before Bayesian network analysis.

### 5. `5_GetMatchedPhenoData.R`
- **Purpose**: Match and merge phenotype data for both datasets.
- **Authors**: Olivia Williamson, Natalie Goulett, Gabriel Odom
- **Details**: Extracts treatment and sample information, transforming treatment data into binary columns representing each treatment group for further analysis.

### 6. `6_BuildBN.R`
- **Purpose**: Build the Bayesian network to analyze gene expression interactions and treatment effects.
- **Authors**: Olivia Williamson, Natalie Goulett, Gabriel Odom
- **Details**: Uses the `bnlearn` package to build a Bayesian network with the Hill Climbing algorithm. Results are visualized and analyzed for significant gene-treatment interactions.

## Datasets

The datasets used in this project:
- **GSE154041**: Contains 74 samples from patients treated with Regorafenib or Lomustine.
- **GSE186332**: Contains 57 samples from a phase-II Selinexor trial with four treatment groups.
  
Both datasets were acquired from [NCBI's GEO database](https://www.ncbi.nlm.nih.gov/geo/).

## Analysis Summary

- **Normalization**: Raw counts were normalized with TMM, followed by transformation to CPM.
- **Merging**: Gene symbols were standardized, and batch effects removed.
- **Bayesian Network Construction**: The Hill Climbing algorithm was used to identify genetic and treatment interactions, highlighting key pathways in GBM progression and treatment response.

## Results and Discussion

The resulting Bayesian network identifies several significant gene-treatment interactions, including the influence of Lomustine on the NGFR gene. Conditional probability tables provide insights into these relationships.

## Dependencies

This project requires the following R packages:
- `BiocManager`, `GEOquery`, `Biobase`, `biomaRt`, `tidyverse`, `hgnc`, `edgeR`, `limma`, `pathwayPCA`, `bnlearn`, `Rgraphviz`

Install packages using:

```R
# Install Bioconductor packages
BiocManager::install(c("GEOquery", "Biobase", "biomaRt", "edgeR", "pathwayPCA", "bnlearn", "Rgraphviz"))

# Install CRAN packages
install.packages(c("tidyverse", "hgnc", "limma"))
