
# 准备环境
rm(list=ls()) 
options(stringsAsFactors = F)
options(timeout = 10000000)
library(Seurat)
library(ggplot2)
library(stringr)
library(tidyverse)
library(dplyr)
library(stringr)
library(data.table)
library(patchwork)
library(ggrepel)
library(CytoTRACE)


load("./8-fibsub/fibsub.Rdata")
load(file = "./0-Rawdata/my_color.Rdata")

## 运行

mat = as.matrix(fibsub@assays$RNA@counts)
phe = as.character(fibsub@meta.data$seurat_clusters)
names(phe) = rownames(fibsub@meta.data)
emb = as.data.frame(fibsub@reductions[["tsne"]]@cell.embeddings)

mat[1:4,1:4]
results <- CytoTRACE(mat = mat)


plotCytoGenes(results, numOfGenes = 10)
plotCytoTRACE(results, gene = "C7", emb = emb)
plotCytoTRACE(results, gene = "MMP11", emb = emb)
plotCytoTRACE(results, gene = "SOX4", emb = emb, outputDir = "./13-CytoTRACE/SOX4.pdf")
plotCytoTRACE(results, phenotype = phe, emb = emb, outputDir = "./13-CytoTRACE/cluster")

plotCytoTRACE(results, phenotype = phe)




