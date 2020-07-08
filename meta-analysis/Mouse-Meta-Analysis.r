library(dplyr)
library(Seurat)
library(patchwork)
library(clustifyr)
library(tidyverse)

mouseAtlas <- readRDS(file.path("~/Reference-Matrix-Generation/atlas/musMusculus/MouseAtlas.rds"))

mouseMetaAnalysis <- CreateSeuratObject(counts = mouseAtlas, project = "Mouse-Meta-Analysis", min.cells = 3, min.features = 200)
mouseMetaAnalysis
gc()

#Normalize Data
#mouseMetaAnalysis <- NormalizeData(mouseMetaAnalysis, normalization.method = "LogNormalize", scale.factor = 10000)

#Preprocessing workflow
mouseMetaAnalysis@assays$RNA@data <- mouseMetaAnalysis@assays$RNA@counts
mouseMetaAnalysis[["percent.mt"]] <- PercentageFeatureSet(mouseMetaAnalysis, pattern = "^MT-")
head(mouseMetaAnalysis@meta.data, 5)
VlnPlot(mouseMetaAnalysis, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot1 <- FeatureScatter(mouseMetaAnalysis, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(mouseMetaAnalysis, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

#Find Variable Features for PCA
mouseMetaAnalysis <- FindVariableFeatures(mouseMetaAnalysis, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(mouseMetaAnalysis), 10)
plot1 <- VariableFeaturePlot(mouseMetaAnalysis)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

#Linear dimension reduction/Run PCA
all.genes <- rownames(mouseMetaAnalysis)
mouseMetaAnalysis <- ScaleData(mouseMetaAnalysis, features = all.genes)
mouseMetaAnalysis <- RunPCA(mouseMetaAnalysis, features = VariableFeatures(object = mouseMetaAnalysis))
print(mouseMetaAnalysis[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(mouseMetaAnalysis, dims = 1:2, reduction = "pca")
DimPlot(mouseMetaAnalysis, reduction = "pca")
DimHeatMap(mouseMetaAnalysis, dims = 1, cells = 500, balanced = TRUE)
DimHeatMap(mouseMetaAnalysis, dims = 1:15, cells = 500, balanced = TRUE)