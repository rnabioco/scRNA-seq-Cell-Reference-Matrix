library(dplyr)
library(Seurat)
library(patchwork)
library(clustifyr)
library(tidyverse)
library(digest)
library(here)

mat_atheroscleroticlesions <- read.table("ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE150nnn/GSE150768/suppl/GSE150768_scRNAseq_rawmatrix.txt.gz", sep = "\t", skipNul = 3)
mat_atheroscleroticlesions <- mat_atheroscleroticlesions %>%
  # as.data.frame() %>%
  column_to_rownames('B11__d1')
#  as.matrix() %>%
#  t()
mat_atheroscleroticlesions[1:5, 1:5]

meta_atheroscleroticlesions <- read_tsv("ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE113nnn/GSE113069/suppl/GSE113069_metadata.txt.gz")
meta_atheroscleroticlesions
ncol(mat_atheroscleroticlesions)

source("~/Reference-Matrix-Generation/R/utils/utils.r")
checkRawCounts(as.matrix(mat_atheroscleroticlesions))

GSE113069Normalized <- NormalizeData(mat_atheroscleroticlesions)

new_ref_matrix <- average_clusters(mat = GSE113069Normalized, metadata = meta_atheroscleroticlesions, if_log = FALSE)
head(new_ref_matrix)
tail(new_ref_matrix)
saveRDS(new_ref_matrix, "GSE150768.rds")
