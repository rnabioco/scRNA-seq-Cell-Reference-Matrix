library(dplyr)
library(Seurat)
library(patchwork)
library(clustifyr)
library(tidyverse)
library(digest)
library(canceR)

GSE129933Filename <- file.choose()
GSE129933Matrix <- readRDS(GSE129933Filename)

GSE137710FilenameMelanoma <- file.choose()
GSE137710Melanoma <- readRDS(GSE137710FilenameMelanoma)

GSE137710FilenameSpleen <- file.choose()
GSE137710Spleen <- readRDS(GSE137710FilenameSpleen)

GSE147405Filename <- file.choose()
GSE147405Matrix <- readRDS(GSE147405Filename)

GSE129933 <- as.data.frame(GSE129933Matrix)
GSE147405 <- as.data.frame(GSE147405Matrix)
GSE137710Melanoma <- as.data.frame(GSE137710Melanoma)
GSE137710Spleen <- as.data.frame(GSE137710Spleen)

rm(GSE129933Matrix)
rm(GSE147405Matrix)

common <- Reduce(intersect, list(rownames(GSE129933), rownames(GSE137710Melanoma), rownames(GSE137710Spleen), rownames(GSE147405)))
GSE129933[common,] # give you common rows in data frame 1  
GSE137710Spleen[common,] # give you common rows in data frame 2
GSE137710Melanoma[common,]
GSE147405[common,]
head(common)

merge(GSE129933, GSE137710Melanoma, GSE137710Spleen, GSE147405, by = "rownames", all = T)
full_join()