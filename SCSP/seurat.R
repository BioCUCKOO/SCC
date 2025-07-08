library(Seurat)
library(dplyr)
library(patchwork)
library(Seurat)
library(data.table)
library(magrittr)
library(DOSE)
library(org.Hs.eg.db)
library(clusterProfiler)
library(tidyverse)
library(tsne)
dev.new()
ori_ct=read.delim("42992_features_clust.txt",sep='\t',row.names = 1,header=T)
data <- t(ori_ct[,-c(1,2,3,4,5)])
maize <- CreateSeuratObject(counts = data,
                            project = "maize")
maize <- NormalizeData(maize)
maize <- FindVariableFeatures(maize, selection.method = "vst", nfeatures = 1500)
dif <- VariableFeatures(maize)

## In addition we scale the data
all.genes <- rownames(maize)
maize <- ScaleData(maize, features = all.genes)
maize <- RunPCA(maize, features = VariableFeatures(object = maize),
                verbose = FALSE)
maize <- FindNeighbors(maize, dims = 1:10, verbose = FALSE)
maize <- FindClusters(maize, resolution =0.8, verbose = FALSE)
maize <- RunUMAP(maize, dims = 1:10, umap.method = "uwot", metric = "cosine")
dplot <- DimPlot(maize, reduction = "umap", group.by = 'seurat_clusters',
                 label = FALSE, pt.size = 2)
dplot
clas <- maize$seurat_clusters
clas <- data.frame(clas)
write.table(clas,"29clust_together_seurat.txt",sep='\t',quote =FALSE)