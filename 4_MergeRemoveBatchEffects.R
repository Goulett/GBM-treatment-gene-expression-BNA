# Title: Merge Data and Remove Batch Effects
# **Olivia Williamson, Natalie Goulett, Gabriel Odom
# Created: 2024-04-22
# Edited: 2024-04-28


# Install and load the following packages:
# BiocManager::install("pathwayPCA")
# install.packages("tidyverse")
# install.packages("limma")
library(tidyverse)
library(pathwayPCA)
library(limma)

GSE154041_Transformed <- read.csv(file = "./Data_Clean/GSE154041_Transformed.csv")
GSE186332_Transformed <- read.csv(file = "./Data_Clean/GSE186332_Transformed.csv")
SharedGenes_Char <- intersect(GSE154041_Transformed$X, GSE186332_Transformed$X)

# # Subset data sects
# SharedGenes_Subset <- inner_join(
#   GSE154041_Transformed,
#   GSE186332_Transformed,
#   by = "X"
# )

# There are 3 duplicated genes: HSPA14, ATXN7, POLR2J3. Let's use R's name
#   repair:
anyDuplicated(GSE154041_Transformed$X)
anyDuplicated(GSE186332_Transformed$X)

###  154041  ###
GSE154041t_df <-
  pathwayPCA::TransposeAssay(GSE154041_Transformed, omeNames = "firstCol") %>%
  as_tibble(.name_repair = "unique")
GSE154041_df <-
  pathwayPCA::TransposeAssay(GSE154041t_df, omeNames = "firstCol")
# Renamed 'X' column to 'gene_symbol'
colnames(GSE154041_df)[1] <- "gene_symbol"

# ONLY RUN ONCE!!
# add dataset name in front of sample id number
# colnames(GSE154041_df)[2:ncol(GSE154041_df)] <-
#   paste0("GSE154041_", colnames(GSE154041_df)[2:ncol(GSE154041_df)])
anyDuplicated(GSE154041_df$gene_symbol)


###  186332  ###
GSE186332t_df <-
  GSE186332_Transformed %>%
  # remove 500-ish genes with missing gene symbols
  filter(!is.na(X)) %>%
  pathwayPCA::TransposeAssay(omeNames = "firstCol") %>%
  # We have three duplicated gene names from the conversion from entrez ID to 
  #   gene symbols. This fixes them.
  as_tibble(.name_repair = "unique")
GSE186332_df <-
  pathwayPCA::TransposeAssay(GSE186332t_df, omeNames = "firstCol")
colnames(GSE186332_df)[1] <- "gene_symbol"
# ONLY RUN ONCE!!
# colnames(GSE186332_df)[2:ncol(GSE186332_df)] <-
#   paste0("GSE186332_", colnames(GSE186332_df)[2:ncol(GSE186332_df)])
anyDuplicated(GSE186332_df$gene_symbol)


###  Inner Join  ###
# we subset data so that each set will have ~13k rows. we have to do inner join
#   so we can use dplyr.
SharedGenes_Subset <- inner_join(
  GSE154041_df,
  GSE186332_df,
  by = "gene_symbol"
)


######  limma  ################################################################
# limma format uses columns for samples. We must create new dataset that is 
#   the inner join of the two sets by column ‘X’.

# We want to separate data into two groups based on trial (not by patient). We 
#   can distinguish between the two using the "_" symbol in the column name.
groups_fct <- colnames(SharedGenes_Subset)[-1] %>%
  # use "_" to distinguish participants from GSE154041
  stringr::str_detect(pattern = "_") %>%
  if_else(true = "GSE154041", false = "GSE186332") %>%
  as_factor()

shared_mat <- as.matrix(SharedGenes_Subset[, -1])
rownames(shared_mat) <- SharedGenes_Subset$gene_symbol

# create new merged dataset without batch effects
Merge_No_BE <- removeBatchEffect(
  x = shared_mat, batch = groups_fct
)

# test_df <- data.frame(
#   y = shared_mat[1,],
#   x = groups_fct
# )
# boxplot(test_df$y ~ test_df$x)

hist(Merge_No_BE)
summary(as.numeric(Merge_No_BE))

# Save
# write.csv(
#   x = Merge_No_BE,
#   file = "Data_Clean/GSE_Data_Joined_Removed_BE_20240426.csv",
#   row.names = TRUE
# )

# Calculate the normal probabilities assoc. w/ threshold z = -1 and z = 1
pnorm(-1)
# 0.1587
qnorm(0.1587)


# Check quantiles to see if it matches the distribution previous:
quantile(as.numeric(Merge_No_BE), probs = c(0.1587, 1 - 0.1587))

# Create function to discretize counts
discretized_genes <- function(normalized_expr_mat,
                              thresh_low = -1,
                              thresh_high = 1){
  # Turn normalized gene expr. values into 3 discrete categories
  # browser()

  # Check for valid threshold values
  if (is.na(thresh_low) | is.na(thresh_high)) {
    stop("NA values in thresholds are not allowed.")
  }

  # Apply cut-offs
  is_low <- normalized_expr_mat < thresh_low
  is_high <- normalized_expr_mat > thresh_high
  is_none <- !(is_low | is_high)

  # Replace z-scores w/ categories
  normalized_expr_mat[is_none] <- 0
  normalized_expr_mat[is_low] <- -1
  normalized_expr_mat[is_high] <- 1

  normalized_expr_mat
}

# Test
# set.seed(1235)
# test <- rnorm(15)
# discretized_genes(test)
# discretized_genes(
#   test, thresh_low = qnorm(0.1587), thresh_high = qnorm(1-0.1587)
# )
#
# # Test on non-normal data
# test2 <- rt(30, df = 3)
# # quantile function returns named vector
# testThresh <- c(
#   low = quantile(test2, 0.1587, names = FALSE),
#   high = quantile(test2, 1 - 0.1587, names = FALSE)
# )
#
# discretized_genes(
#   test2, thresh_low = testThresh["low"], thresh_high = testThresh["high"]
# )

# Test that our error function is working
# discretized_genes(
#   test2, thresh_low = NA, thresh_high = testThresh["high"]
# )

# discretized_genes(
#   test2, thresh_low = testThresh["Low"], thresh_high = testThresh["high"]
# )

# Stopped here
# Batch Effects still - using as example
# discretized_merge_BE <- discretized_genes(
#   shared_mat, thresh_low = qnorm(0.1587), thresh_high = qnorm(1-0.1587)
# )

# After BE removed, the empirical quantiles for 15.87% are +/- .85;
# Should we use these, or should we use +/- 1?
discretized_merge <- discretized_genes(
  Merge_No_BE, thresh_low = -1, thresh_high = 1
  # Merge_No_BE, thresh_low = -0.85, thresh_high = 0.85
)

# We're using compressed data format since most of the data are 0's
# RDS format allows us to compress the size from 24MB to 570KB
saveRDS(discretized_merge, file = "Data_Clean/GSE_discretized_20240426.rds")
