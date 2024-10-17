
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

load(file = "./1-QC/sce.all.filter.Rdata")
load(file = "./0-Rawdata/my_color.Rdata")
color1 = col_vector[1:24]

# 标准化数据
sce.all.filt <- NormalizeData(sce.all.filt, 
                              normalization.method = "LogNormalize",
                              scale.factor = 1e4) 


# 高变异基因的选择
sce.all.filt <- FindVariableFeatures(sce.all.filt)
## Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(sce.all.filt), 10)
## plot variable features with and without labels
plot1 <- VariableFeaturePlot(sce.all.filt)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
ggsave(plot2, filename = "./2-harmony/1-VariableFeatures.pdf", width = 6.25, height = 4.1)

# 归一化数据
## 目的是: 改变每个基因的表达，使细胞的平均表达为0;缩放每个基因的表达，使细胞之间的方差为 1
sce.all.filt <- ScaleData(sce.all.filt)

# PCA分析
sce.all.filt <- RunPCA(sce.all.filt, features = VariableFeatures(object = sce.all.filt))
print(sce.all.filt[["pca"]], dims = 1:5, nfeatures = 5)
p1 = VizDimLoadings(sce.all.filt, dims = 1:4, reduction = "pca", col = 2)
ggsave(p1, filename = "./2-harmony/2-PCA_Vizdim.pdf", width = 8, height = 10)
p2 = DimHeatmap(sce.all.filt, dims = 1:6, cells = 500, balanced = TRUE)
p3 = DimPlot(sce.all.filt, reduction = "pca", group.by = "orig.ident", cols = color1)
p4 = DimPlot(sce.all.filt, reduction = "pca", group.by = "Phase") + scale_color_nejm()
ggsave(p3+p4, filename = "./2-harmony/3-PCA_orig_phase.pdf", width = 10.2, height = 4)

# 确定数据集的"主成分个数"
## NOTE: This process can take a long time for big datasets, comment out for expediency. More
## approximate techniques such as those implemented in ElbowPlot() can be used to reduce
## computation time
sce.all.filt <- JackStraw(sce.all.filt, num.replicate = 100)
sce.all.filt <- ScoreJackStraw(sce.all.filt, dims = 1:20)
JackStrawPlot(sce.all.filt, dims = 1:15)
ElbowPlot(sce.all.filt)

# 整合批次效应
library(harmony)
sce.all.filt <- RunHarmony(sce.all.filt, group.by.vars = "orig.ident")
names(sce.all.filt@reductions)
sce.all.filt <- RunUMAP(sce.all.filt, dims = 1:15, reduction = "harmony")
sce.all.filt <- RunTSNE(sce.all.filt, dims = 1:15, reduction = "harmony")
p1 = DimPlot(sce.all.filt, reduction = "umap", group.by = "orig.ident", label = F, cols = color1)
p2 = DimPlot(sce.all.filt, reduction = "tsne", group.by = "orig.ident", label = F, cols = color1) 
ggsave(p1+p2, filename = "./2-harmony/4-harmony_orig.pdf", width = 11, height = 4.5)
sce.all.filt.harmony = sce.all.filt
save(sce.all.filt.harmony, file = "./3-DoubletFinder/sce.all.filt.harmony.Rdata")


# 细胞聚类
rm(list = ls())
options(stringsAsFactors = F) 
load(file = "./3-DoubletFinder/sce.all.filt.harmony.Rdata")

sce.all.filt.harmony <- FindNeighbors(sce.all.filt.harmony, reduction = "harmony",
                                      dims = 1:15) 

## 设置不同的分辨率，观察分群效果(选择哪一个？)
for (res in c(0.01, 0.05, 0.1, 0.2, 0.3, 0.5, 0.8, 1)) {
  sce.all.filt.harmony = FindClusters(sce.all.filt.harmony, resolution = res, algorithm = 1)
}
colnames(sce.all.filt.harmony@meta.data)
apply(sce.all.filt.harmony@meta.data[,grep("RNA_snn",colnames(sce.all.filt.harmony@meta.data))], 2, table)

FeaturePlot(object = sce.all.filt.harmony, features = "nCount_RNA", reduction = "tsne")


## UMAP查看分群情况
p1_dim = plot_grid(ncol = 4, 
                   DimPlot(sce.all.filt.harmony, reduction = "umap", group.by = "RNA_snn_res.0.01") + ggtitle("resolution_0.01"), 
                   DimPlot(sce.all.filt.harmony, reduction = "umap", group.by = "RNA_snn_res.0.05") + ggtitle("resolution_0.05"), 
                   DimPlot(sce.all.filt.harmony, reduction = "umap", group.by = "RNA_snn_res.0.1") + ggtitle("resolution_0.1"), 
                   DimPlot(sce.all.filt.harmony, reduction = "umap", group.by = "RNA_snn_res.0.2") + ggtitle("resolution_0.2"))
ggsave(plot = p1_dim, filename = "./2-harmony/5-UMAP_resolution_low.pdf", width = 18, height = 4)

p1_dim = plot_grid(ncol = 4, 
                   DimPlot(sce.all.filt.harmony, reduction = "umap", group.by = "RNA_snn_res.0.3") + ggtitle("resolution_0.3"), 
                   DimPlot(sce.all.filt.harmony, reduction = "umap", group.by = "RNA_snn_res.0.5") + ggtitle("resolution_0.5"), 
                   DimPlot(sce.all.filt.harmony, reduction = "umap", group.by = "RNA_snn_res.0.8") + ggtitle("resolution_0.8"), 
                   DimPlot(sce.all.filt.harmony, reduction = "umap", group.by = "RNA_snn_res.1") + ggtitle("resolution_1"))
ggsave(plot = p1_dim, filename = "./2-harmony/6-UMAP_resolution_high.pdf", width = 16, height = 4)


## tSNE查看分群情况
p1_dim = plot_grid(ncol = 4, 
                   DimPlot(sce.all.filt.harmony, reduction = "tsne", group.by = "RNA_snn_res.0.01") + ggtitle("resolution_0.01"), 
                   DimPlot(sce.all.filt.harmony, reduction = "tsne", group.by = "RNA_snn_res.0.05") + ggtitle("resolution_0.05"), 
                   DimPlot(sce.all.filt.harmony, reduction = "tsne", group.by = "RNA_snn_res.0.1") + ggtitle("resolution_0.1"), 
                   DimPlot(sce.all.filt.harmony, reduction = "tsne", group.by = "RNA_snn_res.0.2") + ggtitle("resolution_0.2"))
ggsave(plot = p1_dim, filename = "./2-harmony/tSNE_resolution_low.pdf", width = 18, height = 4)

p1_dim = plot_grid(ncol = 4, 
                   DimPlot(sce.all.filt.harmony, reduction = "tsne", group.by = "RNA_snn_res.0.3") + ggtitle("resolution_0.3"), 
                   DimPlot(sce.all.filt.harmony, reduction = "tsne", group.by = "RNA_snn_res.0.5") + ggtitle("resolution_0.5"), 
                   DimPlot(sce.all.filt.harmony, reduction = "tsne", group.by = "RNA_snn_res.0.8") + ggtitle("resolution_0.8"), 
                   DimPlot(sce.all.filt.harmony, reduction = "tsne", group.by = "RNA_snn_res.1") + ggtitle("resolution_1"))
ggsave(plot = p1_dim, filename = "./2-harmony/tSNE_resolution_high.pdf", width = 20, height = 4)


## 树状图
p2_tree = clustree(sce.all.filt.harmony@meta.data, prefix = "RNA_snn_res.")
ggsave(p2_tree, filename="./2-harmony/7-Tree_diff_resolution.pdf", width = 8.29, height = 8.5)


