# Title: Discretizing the gene expression values
# Author: Olivia Williamson, Natalie Goulett, Gabriel Odom
# Created: 2024-04-18
# Edited: 2024-04-28


# We do not normalize by treatment group. We want to treat all columns as the 
#   same group because this is a preprocessing step. 
# We only use DGEList for raw count data.
# Explanation: https://support.bioconductor.org/p/p133964/#p134019

#install and load edgeR
BiocManager::edgeR
library(edgeR)

# Importing data set files (with gene symbols)
GSE154041 <- read.csv(file = "./Data_Clean/GSE154041_wsymbols_20240410.csv")
GSE186332 <- read.csv(file = "./Data_Clean/GSE186332_wsymbols_20240410.csv")


######  Data Set 1  ############################################################

# GSE186332 - raw counts, converting to normalized counts through TMM method
# Creating DGEList object
test <- DGEList(GSE186332[-2], annotation.columns = "Gene_name")
# **convert something to factors to get CPM to work?
test_TMM <- calcNormFactors(test)
# Convert reads to Counts Per Million (CPM)
test_cpm <- cpm(test_TMM, log = TRUE)
head(test_cpm)
# Check for a relatively normal distribution
hist(test_cpm)
# Transpose matrix to scale, transpose back to long format
test_z <- t(
  scale(
    t(
      test_cpm
    )
  )
)
head(test_z)
hist(test_z)
# There are no hugo gene symbols for 596 of the 25040 genes
rownames(test_z) <- GSE186332$Gene_name
# Cut these to over-, non-, or under-expressed genes


######  Data Set 2  ############################################################

# GSE154041 - TMM Normalized expression values obtained from GEO
# Removing ensemble column
GSE154041_df <- GSE154041[-2]
# Converting to matrix
GSE154041_matrix <- as.matrix(GSE154041_df[-1])
# moving first column to row names - obtain just the values in columns
rownames(GSE154041_matrix) <- GSE154041_df$Gene_name
head(GSE154041_matrix)
# log base 2 (x+1) transformation
GSE154041_log <- log2(GSE154041_matrix + 1)
GSE154041Log2T_mat <- t(
  scale(
    t(
      GSE154041_log
    )
  )
)

hist(GSE154041Log2T_mat)
hist(rnorm(n = 1113990))
# Compared the histogram of our log 2 transformed data set against rnorm
# add more comments/explanation here **lol
# test_z is the transformed version of GSE186332.
# we will read these into the RemoveBatchEffects script.

# Save transformed files
write.csv(
  test_z,
  file = "./Data_Clean/GSE186332_Transformed.csv",
  row.names = TRUE
)
write.csv(
  GSE154041Log2T_mat,
  file = "./Data_Clean/GSE154041_Transformed.csv",
  row.names = TRUE
)

GSE154041_Transformed <- read.csv(
  file = "./Data_Clean/GSE154041_Transformed.csv"
)
GSE186332_Transformed <- read.csv(
  file = "./Data_Clean/GSE186332_Transformed.csv"
)
intersection(GSE154041_Transformed, GSE186332_Transformed)


###### Decompose Normalized Data to remove batch effects ######################
# **are we doing anything with this section or should it be moved? 
#   do you want to keep these links?:
# https://academic.oup.com/bib/article/14/4/469/191565
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3046121/
# https://www.biostars.org/p/266507/
# removeBatchEffects() in Limma::
# because one data set has raw counts, the other TMM normalized (CPM), we can't
#   not using ComBatSeq because one data set is not in raw counts
# after removing BE, data may not come from normal distribution (mean != 0, sd 
#   != 1)
library(limma)

# Merge by gene symbol
intersect(GSE154041Log2T_mat, test_TMM)
removeBatchEffect()


###### Discretizing Values ####################################################

# We need to be able to specify what threshold to use because genomic data is 
#   not normally distributed.
# calculate the normal probabilities assoc. with threshold z = -1 and z = 1
pnorm(-1)
# 0.1587
qnorm(0.1587)

# GSE154041
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
set.seed(1235)
test <- rnorm(15)
discretized_genes(test)
discretized_genes(
  test, thresh_low = qnorm(0.1587), thresh_high = qnorm(1-0.1587)
)

# Test on non-normal data
test2 <- rt(30, df = 3)
# quantile function returns named vector
testThresh <- c(
  low = quantile(test2, 0.1587, names = FALSE),
  high = quantile(test2, 1 - 0.1587, names = FALSE)
)
discretized_genes(
  test2, thresh_low = testThresh["low"], thresh_high = testThresh["high"]
)

# Test that our error function is working
discretized_genes(
  test2, thresh_low = NA, thresh_high = testThresh["high"]
)

discretized_genes(
  test2, thresh_low = testThresh["Low"], thresh_high = testThresh["high"]
)

discrete_counts <- apply(
  standardized.counts,
  MARGIN = 2,
  discretize_genes
)

# Double Check work
discrete_counts[1:5,1:5]
#
standardized.counts[1:5,1:5]
