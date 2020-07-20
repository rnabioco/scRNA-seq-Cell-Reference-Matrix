library(dplyr)
library(Seurat)
library(patchwork)
library(clustifyr)
library(tidyverse)
library(digest)
library(here)

proj_dir <- here()
GSE119228 <- file.path(proj_dir, "Reference-Matrix-Generation", "data", "GSE119228_RAW", "GSM3361719_AB894.txt.gz")

mat_LungDevelopment <- read_tsv(GSE119228)
mat_LungDevelopment <- mat_LungDevelopment %>%
  # as.data.frame() %>%
  column_to_rownames('W171665')
#  as.matrix() %>%
#  t()
mat_LungDevelopment[1:5, 1:5]

meta_LungDevelopment <- read_tsv("ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE119nnn/GSE119228/suppl/GSE119228_metadata.txt.gz")
meta_LungDevelopment <- meta_LungDevelopment %>% column_to_rownames("X1") %>% t() %>% as.data.frame()
meta2 <- meta_LungDevelopment[colnames(mat_LungDevelopment),]
meta_LungDevelopment
ncol(mat_LungDevelopment)

source("~/Reference-Matrix-Generation/R/utils/utils.r")
checkRawCounts(as.matrix(mat_LungDevelopment))

GSE119228Normalized <- NormalizeData(mat_LungDevelopment)

new_ref_matrix <- average_clusters(mat = GSE119228Normalized, metadata = meta_LungDevelopment, cluster_col = "cluster_annotation", if_log = FALSE)
new_ref_matrix_hashed <- average_clusters(mat = mat_LungDevelopment, metadata = meta_LungDevelopment, if_log = FALSE)
head(new_ref_matrix)
tail(new_ref_matrix)
saveRDS(new_ref_matrix, "GSE119228.rds")
