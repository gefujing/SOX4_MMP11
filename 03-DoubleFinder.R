
# 准备环境
rm(list = ls()) 
options(stringsAsFactors = F) 

library(Seurat)
library(ggplot2)
library(clustree)
library(cowplot)
library(dplyr)
library(stringr)
library(data.table)
library(DoubletFinder)

# 准备数据(Seurat标准流程, 标准化、归一化后的数据)
load(file = "./3-DoubletFinder/sce.all.filt.harmony.Rdata")
sce.all.filt.harmony <- FindNeighbors(sce.all.filt.harmony, reduction = "harmony", dims = 1:15) 
sce.all.filt.harmony <- FindClusters(sce.all.filt.harmony, resolution = 0.5)

# 寻找最优pK值
sweep.res.list <- paramSweep_v3(sce.all.filt.harmony, PCs = 1:15, sct = F)
sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)  
bcmvn <- find.pK(sweep.stats)
pK_bcmvn <- bcmvn$pK[which.max(bcmvn$BCmetric)] %>% as.character() %>% as.numeric()

## 排除不能检出的同源doublets，优化期望的doublets数量
DoubletRate = 0.016     # 5000细胞对应的doublets rate是3.9%
homotypic.prop <- modelHomotypic(sce.all.filt.harmony$seurat_clusters)   # 最好提供celltype
nExp_poi <- round(DoubletRate*ncol(sce.all.filt.harmony)) 
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

## 使用确定好的参数鉴定doublets
# 使用nExp = nExp_poi和nExp = nExp_poi.adj,分别进行doublets鉴定，以便后续确定哪些细胞是Doublet-High Confidience
sce.all.filt.harmony <- doubletFinder_v3(sce.all.filt.harmony, PCs = 1:15, 
                                         pN = 0.25, pK = pK_bcmvn,
                                         nExp = nExp_poi, reuse.pANN = F, sct = F)
sce.all.filt.harmony <- doubletFinder_v3(sce.all.filt.harmony, PCs = 1:15, 
                                         pN = 0.25, pK = pK_bcmvn,
                                         nExp = nExp_poi.adj, reuse.pANN = F, sct = F)

## 结果展示
sce.all.filt.harmony@meta.data[,"DF_hi.lo"] <- sce.all.filt.harmony@meta.data$DF.classifications_0.25_0.02_605
sce.all.filt.harmony@meta.data$DF_hi.lo[which(sce.all.filt.harmony@meta.data$DF_hi.lo == "Doublet" & sce.all.filt.harmony@meta.data$DF.classifications_0.25_0.02_539 == "Singlet")] <- "Doublet-Low Confidience"
sce.all.filt.harmony@meta.data$DF_hi.lo[which(sce.all.filt.harmony@meta.data$DF_hi.lo == "Doublet")] <- "Doublet-High Confidience"
table(sce.all.filt.harmony@meta.data$DF_hi.lo)

DimPlot(sce.all.filt.harmony, reduction = "tsne", group.by = "DF_hi.lo", cols =c("gold","red","black"))
ggsave(filename = "./3-DoubletFinder/1-doublet.pdf", width = 5.9, height = 3.8)

sce.har.db = subset(x = sce.all.filt.harmony, DF_hi.lo == "Singlet")
dim(sce.all.filt.harmony)
dim(sce.har.db)
save(sce.har.db, file = "./3-DoubletFinder/sce.har.db.Rdata")
