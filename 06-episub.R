
# 准备环境
rm(list=ls()) 
options(stringsAsFactors = F) 
library(Seurat)
library(ggplot2)
library(stringr)
library(clustree)
library(cowplot)
library(dplyr)
library(stringr)
library(data.table)
library(patchwork)
library(ggsci)
load("./4-celltype/sce.har.db.Rdata")
load(file = "./0-Rawdata/my_color.Rdata")

# 亚群再分析(epi)
table(sce.har.db@meta.data$celltype)
episub = subset(x = sce.har.db, celltype == "Ductal")
table(episub@meta.data$orig.ident)

## 创建Seurat数据
counts = episub@assays$RNA@counts
meta.data = episub@meta.data[, 1:6]
episub = CreateSeuratObject(counts = counts, meta.data = meta.data)
rm(list = c("counts", "meta.data", "sce.har.db"))

## 添加分组信息
table(episub@meta.data$PNI)
table(episub@meta.data$MVI)
table(episub@meta.data$stage)

## 标准化，归一化，PCA
episub <- NormalizeData(episub, normalization.method =  "LogNormalize", scale.factor = 1e4)
episub <- FindVariableFeatures(episub, selection.method = "vst", nfeatures = 2000) 
episub <- ScaleData(episub) 
episub <- RunPCA(object = episub, features = VariableFeatures(episub))
ElbowPlot(episub) 

## tsne
episub <- RunUMAP(episub, dims = 1:20, reduction = "pca")
episub <- RunTSNE(episub, dims = 1:20, reduction = "pca")
p1 = DimPlot(episub, reduction = "umap", group.by = "orig.ident", label = F, cols = col_vector)
ggsave(p1, filename = "./6-episub/1-umap.pdf", width = 5.2, height = 3.9)
p1 = DimPlot(episub, reduction = "tsne", group.by = "orig.ident", label = F, cols = col_vector)
ggsave(p1, filename = "./6-episub/2-tsne.pdf", width = 5.2, height = 3.9)
p1 = DimPlot(episub, reduction = "tsne", group.by = "PNI", label = F, cols = col_vector)
ggsave(p1, filename = "./6-episub/2-tsne.PNI.pdf", width = 5.2, height = 3.9)
p1 = DimPlot(episub, reduction = "tsne", group.by = "MVI", label = F, cols = col_vector)
ggsave(p1, filename = "./6-episub/2-tsne.MVI.pdf", width = 5.2, height = 3.9)
p1 = DimPlot(episub, reduction = "tsne", group.by = "stage", label = F, cols = col_vector)
ggsave(p1, filename = "./6-episub/2-tsne.stage.pdf", width = 5.2, height = 3.9)


## 细胞聚类
episub <- FindNeighbors(episub, reduction = "pca", dims = 1:20) 
episub = FindClusters(episub, resolution = 0.01, algorithm = 1)
colnames(episub@meta.data)
table(episub@meta.data$seurat_clusters)

# episub = subset(episub, RNA_snn_res.0.01 %in% c(0,1,2))

## tsne查看分群情况
p1_dim = DimPlot(episub, reduction = "tsne", group.by = "seurat_clusters", label = T, cols = col_vector)
ggsave(plot = p1_dim, filename = "./6-episub/3-epi.cluster.umap.pdf", width = 3.3, height = 2.8)


## 组别比例
## 绘制堆叠条形图(PNI)
cell.prop = as.data.frame(prop.table(table(episub$seurat_clusters, episub$PNI), 2))
colnames(cell.prop) = c("Celltype", "Group", "Proportion")
cell.prop$proportion2 = round(cell.prop$Proportion, 2)

th = theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0.5)) 
p1 = ggplot(cell.prop, aes(Group, Proportion, fill = Celltype)) +
  geom_bar(stat = "identity", position = "fill") +
  ggtitle("Celltypes of Group") +
  theme_classic() +
  guides(fill = guide_legend(title = NULL)) + th + scale_fill_manual(values = col_vector)

## 绘制堆叠条形图(MVI)
cell.prop = as.data.frame(prop.table(table(episub$seurat_clusters, episub$MVI), 2))
colnames(cell.prop) = c("Celltype", "Group", "Proportion")
cell.prop$proportion2 = round(cell.prop$Proportion, 2)

th = theme(axis.text.x = element_text(angle = 90, 
                                      vjust = 0.5, hjust=0.5)) 
p2 = ggplot(cell.prop, aes(Group, Proportion, fill = Celltype)) +
  geom_bar(stat = "identity", position = "fill") +
  ggtitle("Celltypes of Group") +
  theme_classic() +
  guides(fill = guide_legend(title = NULL)) + th + scale_fill_manual(values = col_vector)

## 绘制堆叠条形图(stage)
cell.prop = as.data.frame(prop.table(table(episub$seurat_clusters, episub$stage), 2))
colnames(cell.prop) = c("Celltype", "Group", "Proportion")
cell.prop$proportion2 = round(cell.prop$Proportion, 2)

th = theme(axis.text.x = element_text(angle = 90, 
                                      vjust = 0.5, hjust=0.5)) 
p3 = ggplot(cell.prop, aes(Group, Proportion, fill = Celltype)) +
  geom_bar(stat = "identity", position = "fill") +
  ggtitle("Celltypes of Group") +
  theme_classic() +
  guides(fill = guide_legend(title = NULL)) + th + scale_fill_manual(values = col_vector)

pall = p1 + p2 + p3 + plot_layout(guides = 'collect')
ggsave(pall, filename = "./6-episub/4-cell.proportion.pdf", width = 5.8, height = 3.4)

save(episub, file = "./6-episub/episub.Rdata")


# 设置多线程
rm(list=ls()) 
load(file = "./6-episub/episub.Rdata")

library(future) # check the current active plan
plan()
plan("multiprocess", workers = 10)
plan()

# 计算marker基因
Idents(episub) = episub@meta.data$seurat_clusters
table(Idents(episub))
sce.markers <- FindAllMarkers(object = episub, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25)
DT::datatable(sce.markers)
write.csv(sce.markers, file = "./6-episub/episub.markers.csv")


