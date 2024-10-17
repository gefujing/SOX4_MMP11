
# 设置环境
rm(list = ls()) 
options(stringsAsFactors = F) 

library(Seurat)
library(ggplot2)
library(clustree)
library(cowplot)
library(dplyr)
library(stringr)
library(data.table)
library(pheatmap)
library(ggheatmap)
library(RColorBrewer)
library(ggthemes)
library(patchwork)
library(scRNAtoolVis)

# 设置数据
load(file = "./4-celltype/sce.har.db.Rdata")
load(file = "./0-Rawdata/my_color.Rdata")
table(sce.har.db$celltype)
table(sce.har.db$DF_hi.lo)


# 单细胞特征
## 细胞标志基因-文献
genes_to_check = c('KRT19','PRSS1', 'CHGB', 'CDH5', 'LUM', 'RGS5', 'AIF1', 'CD3D', 'MS4A1', 'MZB1')
p <- DotPlot(sce.har.db, features = unique(genes_to_check), group.by = "celltype", assay='RNA')
ggsave(filename = "./5-celltype-feature/2-celltype.markers.pdf", width = 9, height = 3.3)

p <- VlnPlot(sce.har.db, features = unique(genes_to_check), group.by = "celltype", assay='RNA', pt.size = 0, ncol = 5, cols = col_vector)
ggsave(filename = "./5-celltype-feature/2-celltype.markers2.pdf", width = 13, height = 5.5)


## 组别比例
## 绘制堆叠条形图(PNI)
cell.prop = as.data.frame(prop.table(table(sce.har.db$celltype, sce.har.db$PNI), 2))
colnames(cell.prop) = c("Celltype", "Group", "Proportion")
cell.prop$proportion2 = round(cell.prop$Proportion, 2)

th = theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0.5)) 
p1 = ggplot(cell.prop, aes(Group, Proportion, fill = Celltype)) +
  geom_bar(stat = "identity", position = "fill") +
  ggtitle("Celltypes of Group") +
  theme_classic() +
  guides(fill = guide_legend(title = NULL)) + th + scale_fill_aaas()

## 绘制堆叠条形图(MVI)
cell.prop = as.data.frame(prop.table(table(sce.har.db$celltype, sce.har.db$MVI), 2))
colnames(cell.prop) = c("Celltype", "Group", "Proportion")
cell.prop$proportion2 = round(cell.prop$Proportion, 2)

th = theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0.5)) 
p2 = ggplot(cell.prop, aes(Group, Proportion, fill = Celltype)) +
  geom_bar(stat = "identity", position = "fill") +
  ggtitle("Celltypes of Group") +
  theme_classic() +
  guides(fill = guide_legend(title = NULL)) + th + scale_fill_aaas()

## 绘制堆叠条形图(STAGE)
cell.prop = as.data.frame(prop.table(table(sce.har.db$celltype, sce.har.db$stage), 2))
colnames(cell.prop) = c("Celltype", "Group", "Proportion")
cell.prop$proportion2 = round(cell.prop$Proportion, 2)

th = theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0.5)) 
p3 = ggplot(cell.prop, aes(Group, Proportion, fill = Celltype)) +
  geom_bar(stat = "identity", position = "fill") +
  ggtitle("Celltypes of Group") +
  theme_classic() +
  guides(fill = guide_legend(title = NULL)) + th + scale_fill_aaas()

## 拼图
pall = p1 + p2 + p3 + plot_layout(guides = 'collect')
ggsave(pall, filename = "./5-celltype-feature/3-cell.proportion.pdf", width = 5.8, height = 3.4)


## 细胞标志基因-find
rm(list = ls()) 
load(file = "./4-celltype/sce.har.db.Rdata")

### 设置多线程
library(future) # check the current active plan
plan()
plan("multiprocess", workers = 10)
plan()

# 计算marker基因
sce.har.db
Idents(sce.har.db) = sce.har.db@meta.data$celltype
table(Idents(sce.har.db))
table(sce.har.db@meta.data$orig.ident)
sce.markers <- FindAllMarkers(object = sce.har.db, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25)
DT::datatable(sce.markers)
write.csv(sce.markers, file = "./5-celltype-feature/cell.marker.csv")

# marker基因可视化
top3 <- sce.markers %>% group_by(cluster) %>% top_n(3, avg_log2FC)
DoHeatmap(sce.har.db, top3$gene[-27], group.by = 'celltype', size = 5)
ggsave(filename = "./5-celltype-feature/4-find.markers.top3.heatmap.pdf")

th = theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0.5)) 
p <- DotPlot(sce.har.db, features = unique(top3$gene), group.by = 'celltype', assay='RNA') + coord_flip() + th
p
ggsave(filename = "./5-celltype-feature/4-find.markers_top3_dotplot.pdf",  width = 5.5, height = 6.7)


# 新式热图1
load(file = "./4-celltype/sce.har.db.Rdata")
table(sce.har.db$celltype)
table(sce.har.db$PNI)
table(sce.har.db$MVI)
table(sce.har.db$stage)
sce.markers = read.csv(file = "./5-celltype-feature/cell.marker.csv", header = T)
rownames(sce.markers) = sce.markers$X
sce.markers = sce.markers[,-1]

deg.sig = sce.markers %>% filter(pct.1 > 0.7 & p_val_adj < 0.001) %>% arrange(cluster, p_val_adj)
deg.top3 = deg.sig %>% group_by(cluster)  %>% top_n(3, avg_log2FC) %>% top_n(3, p_val_adj) %>% arrange(cluster, p_val_adj)

avgData = sce.har.db@assays$RNA@data[deg.top3$gene, ] %>% 
  apply(1, function(x){
    tapply(x, sce.har.db$celltype, mean)
  }) %>% t

phData = MinMax(scale(avgData),-2,2)
# rownames(phData) = 1:nrow(phData)

rownames(deg.top3) = deg.top3$gene
annotation_row =  data.frame(cluster = deg.top3$cluster)
rownames(annotation_row) = rownames(deg.top3)
phData = phData[,c(2,9,1,3,4,5,7,10,6,8)]

pheatmap(phData,
         color = colorRampPalette(c("darkblue", "white", "red3"))(99),
         scale = "row",
         cluster_cols = F,
         cluster_rows = F,
         annotation_row = annotation_row
         )

# 新式热图2
load(file = "./4-celltype/sce.har.db.Rdata")
table(sce.har.db$celltype)

sce.markers = read.csv(file = "./5-celltype-feature/cell.marker.csv", header = T)
rownames(sce.markers) = sce.markers$X
sce.markers = sce.markers[,-1]

Idents(sce.har.db) = "celltype"
sce.har.db = subset(sce.har.db, downsample = 100)

top5 = sce.markers %>% group_by(cluster) %>% top_n(n=5, wt = avg_log2FC)

colanno = sce.har.db@meta.data[,c("celltype", "PNI", "MVI", "stage")]
colanno = colanno %>% arrange(celltype)
colanno$celltype = factor(colanno$celltype, levels = unique(colanno$celltype))

rowanno = top5
# rowanno = rowanno %>% arrange(cluster)

dat.exp = sce.har.db@assays$RNA@counts[rowanno$gene, rownames(colanno)]
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


# 新式气泡图
top3 = sce.markers %>% group_by(cluster) %>% top_n(n=3, wt = avg_log2FC)
dat.plot = as.matrix(sce.har.db@assays$RNA@data[top3$gene,])
dat.plot = t(dat.plot)
dat.plot = as.data.frame(scale(dat.plot))
dat.plot$barcode = rownames(dat.plot)

anno = sce.har.db@meta.data
anno = anno[,c("orig.ident", "celltype")]
anno$barcode = rownames(anno)
dat.plot = inner_join(x = dat.plot, y = anno, by = "barcode")
dat.plot$barcode = NULL

celltype_v=c()
gene_v=c()
mean_v=c()
ratio_v=c()
for (i in unique(dat.plot$celltype)) {
  dat.plot_small = dat.plot %>% filter(celltype == i)
  for (j in top3$gene) {
    exp_mean = mean(dat.plot_small[,j])
    exp_ratio = sum(dat.plot_small[,j] > min(dat.plot_small[,j])) / length(dat.plot_small[,j])
    celltype_v = append(celltype_v,i)  ##Add elements to a vector.
    gene_v = append(gene_v,j)
    mean_v = append(mean_v,exp_mean)
    ratio_v = append(ratio_v,exp_ratio)
  }
}

plotdf = data.frame(celltype = celltype_v,
                    gene = gene_v,
                    exp = mean_v,
                    ratio = ratio_v)

plotdf$celltype = factor(plotdf$celltype, levels = c('Ductal', 'Acinar', 'Endocrine', 'Endothelial', 
                                                     'Fibroblast', 'Stellate', 'Macrophage', 
                                                     'T cell', 'B cell', 'Plasma cell'))
plotdf$gene = factor(plotdf$gene, levels = rev(as.character(top5$gene)))
plotdf$exp = ifelse(plotdf$exp > 3, 3, plotdf$exp)

library(RColorBrewer)
plotdf %>% 
  ggplot(aes(x = celltype, y=gene, size = ratio, color = exp)) + 
  geom_point() +
  scale_x_discrete("") +
  scale_y_discrete("") +
  scale_color_gradientn(colours = rev(c("#FFD92F","#FEE391",brewer.pal(11, "Spectral")[7:11]))) +
  scale_size_continuous(limits = c(0,1)) +
  theme_bw() +
  theme(axis.text.x.bottom = element_text(hjust = 1, vjust = 1, angle = 45))
ggsave(filename = "./5-celltype-feature/6-celltyp-dotplot.pdf", width = 4.5, height = 6.5)


# 新式特征图
genes_to_check = c('KRT19','PRSS1', 'CHGB', 'CDH5', 'LUM',
                   'RGS5', 'AIF1', 'CD3D', 'MS4A1')

FeaturePlot(object = sce.har.db, 
            features = "SEC14L2", 
            reduction = "tsne", 
            ncol = 1,
            cols = c("grey", "red"))
ggsave(filename = "./5-celltype-feature/6-celltype.feature.pdf", width = 14, height = 5)

DimPlot(sce.har.db, group.by = "celltype", reduction = "tsne")


# 火山图
load(file = "./4-celltype/sce.har.db.Rdata")
load(file = "./0-Rawdata/my_color.Rdata")
table(sce.har.db$celltype)
Idents(sce.har.db) = "celltype"

sce.markers = read.csv(file = "./5-celltype-feature/cell.marker.csv", header = T)
rownames(sce.markers) = sce.markers$X
sce.markers = sce.markers[,-1]
top5 = sce.markers %>% group_by(cluster) %>% top_n(n=5, wt = avg_log2FC)

p = AverageHeatmap(object = sce.har.db, markerGene = top5$gene, annoCol = T, myanCol = col_vector[1:10]) 



# 新式特征图

FeaturePlot(object = sce.har.db, 
            features = "SOX4", 
            reduction = "tsne", 
            ncol = 1,
            cols = c("grey", "red"))
ggsave(filename = "./5-celltype-feature/6-celltype.SOX4.feature.pdf", width = 3.8, height = 3.3)






