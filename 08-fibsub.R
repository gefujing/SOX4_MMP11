
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
library(RColorBrewer)
library(harmony)
load("./4-celltype/sce.har.db.Rdata")
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
#处理后有73种差异还比较明显的颜色，基本够用
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals))) 


# 亚群再分析(fib)
table(sce.har.db@meta.data$celltype)
fibsub = subset(x = sce.har.db, celltype == "Fibroblast")

## 创建Seurat数据
counts = fibsub@assays$RNA@counts
meta.data = fibsub@meta.data[, 1:6]
fibsub = CreateSeuratObject(counts = counts, meta.data = meta.data)
rm(list = c("counts", "meta.data", "sce.har.db"))

## 标准化，归一化，PCA
fibsub <- NormalizeData(fibsub, normalization.method =  "LogNormalize", scale.factor = 1e4)
fibsub <- FindVariableFeatures(fibsub, selection.method = "vst", nfeatures = 2000) 
fibsub <- ScaleData(fibsub) 
fibsub <- RunPCA(object = fibsub, features = VariableFeatures(fibsub))
DimHeatmap(fibsub, dims = 1:12, cells = 100, balanced = TRUE)
ElbowPlot(fibsub) 

## harmony
fibsub <- RunHarmony(fibsub, group.by.vars = "orig.ident")
names(fibsub@reductions)
fibsub <- RunUMAP(fibsub, dims = 1:20, reduction = "harmony")
fibsub <- RunTSNE(fibsub, dims = 1:20, reduction = "harmony")
p1 = DimPlot(fibsub, reduction = "umap", group.by = "orig.ident", label = F, cols = col_vector)
ggsave(p1, filename = "./8-fibsub/1-harmony.umap.pdf", width = 5.2, height = 3.9)
p1 = DimPlot(fibsub, reduction = "tsne", group.by = "orig.ident", label = F, cols = col_vector)
ggsave(p1, filename = "./8-fibsub/2-harmony.tsne.pdf", width = 5.2, height = 3.9)


## 细胞聚类
fibsub = FindNeighbors(fibsub, reduction = "harmony", dims = 1:20) 
fibsub = FindClusters(fibsub, resolution = 0.1, algorithm = 1)
colnames(fibsub@meta.data)
table(fibsub@meta.data$seurat_clusters)

fibsub = subset(fibsub, seurat_clusters %in% c(0,1))

## tsne查看分群情况
p1_dim = DimPlot(fibsub, reduction = "tsne", group.by = "seurat_clusters", label = T, cols = col_vector)
ggsave(plot = p1_dim, filename = "./8-fibsub/3-epi.cluster.umap.pdf", width = 4.5, height = 4)


## 组别比例
## 绘制堆叠条形图-PNI
cell.prop = as.data.frame(prop.table(table(fibsub$seurat_clusters, fibsub$PNI), 2))
colnames(cell.prop) = c("Celltype", "Group", "Proportion")
cell.prop$proportion2 = round(cell.prop$Proportion, 2)

th = theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0.5)) 
p1 = ggplot(cell.prop, aes(Group, Proportion, fill = Celltype)) +
  geom_bar(stat = "identity", position = "fill") +
  ggtitle("Celltypes of Group") +
  theme_classic() +
  guides(fill = guide_legend(title = NULL)) + scale_fill_nejm()


## 绘制堆叠条形图-MVI
cell.prop = as.data.frame(prop.table(table(fibsub$seurat_clusters, fibsub$MVI), 2))
colnames(cell.prop) = c("Celltype", "Group", "Proportion")
cell.prop$proportion2 = round(cell.prop$Proportion, 2)

p2 = ggplot(cell.prop, aes(Group, Proportion, fill = Celltype)) +
  geom_bar(stat = "identity", position = "fill") +
  ggtitle("Celltypes of Group") +
  theme_classic() +
  guides(fill = guide_legend(title = NULL)) + scale_fill_nejm()


## 绘制堆叠条形图-stage
cell.prop = as.data.frame(prop.table(table(fibsub$seurat_clusters, fibsub$stage), 2))
colnames(cell.prop) = c("Celltype", "Group", "Proportion")
cell.prop$proportion2 = round(cell.prop$Proportion, 2)

p3 = ggplot(cell.prop, aes(Group, Proportion, fill = Celltype)) +
  geom_bar(stat = "identity", position = "fill") +
  ggtitle("Celltypes of Group") +
  theme_classic() +
  guides(fill = guide_legend(title = NULL)) + scale_fill_nejm()

pall = p1 + p2 + p3 + plot_layout(guides = 'collect')
ggsave(pall, filename = "./8-fibsub/4-cell.proportion.pdf", width = 5.8, height = 3.4)

save(fibsub, file = "./8-fibsub/fibsub.Rdata")


# 设置多线程
library(future) # check the current active plan
plan()
plan("multiprocess", workers = 10)
plan()

# 计算marker基因
Idents(fibsub) = fibsub@meta.data$seurat_clusters
table(Idents(fibsub))
table(fibsub@meta.data$orig.ident)
sce.markers <- FindAllMarkers(object = fibsub, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25)
DT::datatable(sce.markers)
write.csv(sce.markers, file = "./8-fibsub/fibsub.markers.csv")

# 热图
fibsub = subset(fibsub, downsample = 500)
top5 = sce.markers %>% group_by(cluster) %>% top_n(n=5, wt = avg_log2FC)

colanno = fibsub@meta.data[,c( "PNI", "MVI", "stage", "seurat_clusters")]
colanno = colanno %>% arrange(seurat_clusters)
colanno$seurat_clusters = factor(colanno$seurat_clusters, levels = unique(colanno$seurat_clusters))

rowanno = top5
# rowanno = rowanno %>% arrange(cluster)

dat.exp = fibsub@assays$RNA@counts[rowanno$gene, rownames(colanno)]
dat.exp = scale(dat.exp)
dat.exp[dat.exp>=2] = 2
dat.exp[dat.exp<=(-2)] = -2

pheatmap(dat.exp,
         color = colorRampPalette(c("darkblue", "white", "red3"))(99),
         cluster_rows = F,
         cluster_cols = F,
         show_colnames = F,
         annotation_col = colanno,
         gaps_row = as.numeric(cumsum(table(rowanno$cluster))[-10]),
         gaps_col = as.numeric(cumsum(table(colanno$celltype))[-10])
)

# 气泡图
top5 = sce.markers %>% group_by(cluster) %>% top_n(n=5, wt = avg_log2FC)
dat.plot = as.matrix(fibsub@assays$RNA@data[top5$gene,])
dat.plot = t(dat.plot)
dat.plot = as.data.frame(scale(dat.plot))
dat.plot$barcode = rownames(dat.plot)

anno = fibsub@meta.data
anno = anno[,c("orig.ident", "seurat_clusters")]
anno$barcode = rownames(anno)
dat.plot = inner_join(x = dat.plot, y = anno, by = "barcode")
dat.plot$barcode = NULL

celltype_v=c()
gene_v=c()
mean_v=c()
ratio_v=c()
for (i in unique(dat.plot$seurat_clusters)) {
  dat.plot_small = dat.plot %>% filter(seurat_clusters == i)
  for (j in top5$gene) {
    exp_mean = mean(dat.plot_small[,j])
    exp_ratio = sum(dat.plot_small[,j] > min(dat.plot_small[,j])) / length(dat.plot_small[,j])
    celltype_v = append(celltype_v,i)  ##Add elements to a vector.
    gene_v = append(gene_v,j)
    mean_v = append(mean_v,exp_mean)
    ratio_v = append(ratio_v,exp_ratio)
  }
}

plotdf = data.frame(seurat_clusters = celltype_v,
                    gene = gene_v,
                    exp = mean_v,
                    ratio = ratio_v)


plotdf$gene = factor(plotdf$gene, levels = rev(as.character(top5$gene)))
plotdf$exp = ifelse(plotdf$exp > 3, 3, plotdf$exp)

library(RColorBrewer)
plotdf %>% 
  ggplot(aes(x = seurat_clusters, y=gene, size = ratio, color = exp)) + 
  geom_point() +
  scale_x_discrete("") +
  scale_y_discrete("") +
  scale_color_gradientn(colours = rev(c("#FFD92F","#FEE391",brewer.pal(11, "Spectral")[7:11]))) +
  scale_size_continuous(limits = c(0,1)) +
  theme_bw() +
  theme(axis.text.x.bottom = element_text(hjust = 1, vjust = 1, angle = 45))
ggsave(filename = "./8-fibsub/6-fibsub-markers.pdf", width = 2.4, height = 3)

