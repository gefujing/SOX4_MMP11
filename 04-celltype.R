
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
library(ggheatmap)
library(ggsci)

# 设置数据
load(file = "./3-DoubletFinder/sce.har.db.Rdata")
load(file = "./0-Rawdata/my_color.Rdata")
color2 = col_vector[1:22]
color1 = col_vector[1:24]


# 检测分群相关性
table(sce.har.db$seurat_clusters)
av <- AverageExpression(sce.har.db, group.by = "seurat_clusters", assays = "RNA")
av <- av[[1]]
cg <- names(tail(sort(apply(av, 1, sd)),1000))
View(av[cg,])
View(cor(av[cg,],method = 'spearman'))
ggheatmap(cor(av[cg,],method = 'spearman'),
          cluster_rows = T,
          cluster_cols = T,
          border = "black",
          color = c("blue", "white", "red"))
ggsave(filename = "./4-celltype/1-cluster.cor.heatmap.pdf", width = 5, height = 4)

DimPlot(object = sce.har.db, group.by = "seurat_clusters", reduction = "tsne", cols = color2) 
ggsave(file = "./4-celltype/2-tsne.cluster.pdf", width = 5, height = 4)

DimPlot(object = sce.har.db, group.by = "orig.ident", reduction = "tsne", cols = color1) 
ggsave(file = "./4-celltype/2-tsne.patient.pdf", width = 5, height = 4)

DimPlot(object = sce.har.db, group.by = "PNI", reduction = "tsne") + scale_color_nejm()
ggsave(file = "./4-celltype/2-tsne.PNI.pdf", width = 5, height = 4)

DimPlot(object = sce.har.db, group.by = "MVI", reduction = "tsne") + scale_color_nejm()
ggsave(file = "./4-celltype/2-tsne.MVI.pdf", width = 5, height = 4)

DimPlot(object = sce.har.db, group.by = "stage", reduction = "tsne") + scale_color_nejm()
ggsave(file = "./4-celltype/2-tsne.stage.pdf", width = 5, height = 4)


# 检查总体基因情况
genes_to_check = c('PTPRC', # immune cell
                   'CD3D', 'CD3E', 'CD4','CD8A', # T Cells 
                   'CD19', 'CD79A', 'MS4A1', # B cells
                   'IGHG1', 'MZB1', 'SDC1', # Plasma cells
                   'CD68', 'CD163', 'CD14', 'C1QA', 'C1QB', 'ITGAM', 'AIF1',# macrophages
                   'TPSAB1', 'TPSB2', # mast cells,
                   'RGS5', 'CD73', 'CD105', 'CD44', # perivascular cellhttp://biotrainee.vip:12133/graphics/plot_zoom_png?width=561&height=787
                   'CD14', 'S100A9', 'S100A8', 'MMP19', # monocyte
                   'FCGR3A', 'FGFBP2', 'CX3CR1', 'KLRB1', 'NCR1', # NK cells
                   'LAMP3', 'IDO1','IDO2',## DC3 
                   'CD1E','CD1C', # DC2
                   'FGF7','MME', 'ACTA2', ## human Fibroblasts 
                   'DCN', 'LUM', 'GSN' , ## mouse PDAC Fibroblasts 
                   'MKI67' , 'TOP2A', 
                   'PECAM1', 'VWF',  'CDH5', ## Endothelial cells
                   'AMY1', 'AMY2A2', 'PRSS1',  ## Acinar cells
                   'EPCAM' , 'KRT19', 'KRT7', 'PROM1', 'ALDH1A1', 'CD24', # epithelial or tumor
                   'CHGB' ## Endocrine cells
)

genes_to_check
p <- DotPlot(sce.har.db, features = unique(genes_to_check), group.by = "seurat_clusters",
             assay='RNA') + coord_flip()
p
ggsave(filename = "./4-celltype/3-cluster.markers.pdf", width = 8, height = 10)


# 检查总体基因情况
genes_to_check = c('PTPRC', # immune cell
                   'CD3D', 'CD3E', # T Cells 
                   'CD79A', 'MS4A1', # B cells
                   'MZB1', 'SDC1', # Plasma cells
                   'CD68', 'CD163', 'CD14', 'C1QA', 'C1QB', 'ITGAM', 'AIF1',# macrophages
                   'RGS5', 'CD44', # perivascular
                   'CD14', 'S100A9', 'S100A8', # monocyte
                   'FCGR3A', # NK cells
                   'CD1E','CD1C', # DC2
                   'FGF7','MME', 'ACTA2', ## human Fibroblasts 
                   'MKI67' , 'TOP2A', 
                   'VWF',  'CDH5', ## Endothelial cells
                   'PRSS1',  ## Acinar cells
                   'EPCAM' , 'KRT19', 'KRT7', 'PROM1', 'ALDH1A1', # epithelial or tumor
                   'CHGB' ## Endocrine cells
)

genes_to_check
p <- DotPlot(sce.har.db, features = unique(genes_to_check), group.by = "seurat_clusters",
             assay='RNA') + coord_flip()
p
ggsave(filename = "./4-celltype/3-cluster.markers2.pdf", width = 8, height = 8)

# 检查总体基因情况
genes_to_check = c('KRT19','PRSS1', 'CHGB', 'CDH5', 'LUM',
                   'RGS5', 'AIF1', 'CD3D', 'MS4A1')

genes_to_check
p <- DotPlot(sce.har.db, features = unique(genes_to_check), group.by = "seurat_clusters",
             assay='RNA') + coord_flip()
p
ggsave(filename = "./4-celltype/3-cluster.markers3.pdf", width = 9.5, height = 3.5)

# 细胞注释
celltype = data.frame(ClusterID = 0:21,
                      celltype = 0:21) 
## 定义细胞亚群
celltype[celltype$ClusterID %in% c(0,5,9,15,17), 2] = 'Ductal'
celltype[celltype$ClusterID %in% c(10), 2] = 'Acinar' 
celltype[celltype$ClusterID %in% c(12), 2] = 'Endocrine'
celltype[celltype$ClusterID %in% c(1,20,21), 2] = 'Endothelial'
celltype[celltype$ClusterID %in% c(4,6,11,16), 2] = 'Fibroblast'
celltype[celltype$ClusterID %in% c(2,18,19), 2] = 'Stellate' 
celltype[celltype$ClusterID %in% c(3,14), 2] = 'Macrophage'
celltype[celltype$ClusterID %in% c(7), 2] = 'T cell'
celltype[celltype$ClusterID %in% c(8), 2] = 'B cell'  
celltype[celltype$ClusterID %in% c(13), 2] = 'Plasma cell'  

## 写入细胞亚群
table(celltype$celltype)
sce.har.db@meta.data$celltype = "NA"

for(i in 1:nrow(celltype)){
  sce.har.db@meta.data[which(sce.har.db@meta.data$seurat_clusters == celltype$ClusterID[i]), 'celltype'] <- celltype$celltype[i]}

table(sce.har.db@meta.data$celltype)
sce.har.db@meta.data$celltype = factor(sce.har.db@meta.data$celltype,
                                         levels = c('Ductal', 'Acinar', 'Endocrine', 'Endothelial', 
                                                    'Fibroblast', 'Stellate', 'Macrophage', 
                                                    'T cell', 'B cell', 'Plasma cell'))
save(sce.har.db, file = "./4-celltype/sce.har.db.Rdata")

# 查看细胞亚群
DimPlot(sce.har.db, reduction = "tsne", group.by = "celltype", label = T, cols = col_vector) 
ggsave(filename = "./4-celltype/4-tsne.celltype.pdf", width = 5, height = 4)

DimPlot(sce.har.db, reduction = "tsne", group.by = "celltype", split.by = "MVI", label = T, cols = col_vector) 
ggsave(filename = "./4-celltype/4-tsne.celltype.MVI.pdf", width = 8, height = 4)

DimPlot(sce.har.db, reduction = "tsne", group.by = "celltype", split.by = "PNI", label = T, cols = col_vector) 
ggsave(filename = "./4-celltype/4-tsne.celltype.PNI.pdf", width = 8, height = 4)

DimPlot(sce.har.db, reduction = "tsne", group.by = "celltype", split.by = "stage", label = T, cols = col_vector) 
ggsave(filename = "./4-celltype/4-tsne.celltype.stage.pdf", width = 8, height = 4)
