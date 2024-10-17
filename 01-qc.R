
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

# 01导入数据
## 读取矩阵，创建Seurat对象
rawcount = fread(file = "./0-Rawdata/count-matrix.txt", data.table = F)
load(file = "./0-Rawdata/my_color.Rdata")
rawcount[1:4, 1:4]
rownames(rawcount) = rawcount[ ,1]
rawcount = rawcount[ ,-1]
sce.all = CreateSeuratObject(counts = rawcount)

## 检查表达矩阵
as.data.frame(sce.all@assays$RNA@counts[1:10, 1:2])
head(sce.all@meta.data, 20)
tail(sce.all@meta.data, 20)
table(sce.all@meta.data$orig.ident) 

## 取子集
Tumor = c("T2", "T3", "T4", "T5", "T6", "T7", "T9", "T10", "T11", "T13", "T14", "T15", 
          "T16", "T17", "T18", "T19", "T20", "T21", "T22", "T23", "T24")
sce.all = subset(x = sce.all, idents = Tumor)
sce.all@meta.data$orig.ident = factor(sce.all@meta.data$orig.ident, 
                                      levels = c("T2", "T3", "T4", "T5", "T6", "T7", "T9", "T10", "T11", "T13", "T14", 
                                                 "T15", "T16", "T17", "T18", "T19", "T20", "T21", "T22", "T23", "T24"))

## 添加分组信息
sjqx = c("T2", "T9", "T10", "T11", "T13", "T14", "T15", "T16", "T17", "T18", "T19", "T20", "T21", "T22", "T23", "T24")
sce.all@meta.data$PNI = ifelse(sce.all@meta.data$orig.ident %in% sjqx, "PNI", "None_PNI")
sce.all@meta.data$PNI = factor(sce.all@meta.data$PNI, levels = c("None_PNI", "PNI"))

xgqx = c("T9","T10","T11","T14","T16", "T18", "T19", "T20", "T21", "T23")
sce.all@meta.data$MVI = ifelse(sce.all@meta.data$orig.ident %in% xgqx, "MVI", "None_MVI")
sce.all@meta.data$MVI = factor(sce.all@meta.data$MVI, levels = c("None_MVI", "MVI"))

stage1 = c("T3","T5","T10","T17","T18","T19","T21","T22","T24")
sce.all@meta.data$stage = ifelse(sce.all@meta.data$orig.ident %in% stage1, "Stage_I", "Stage_II")
sce.all@meta.data$stage = factor(sce.all@meta.data$stage, levels = c("Stage_I", "Stage_II"))

table(sce.all@meta.data$PNI)
table(sce.all@meta.data$MVI)
table(sce.all@meta.data$stage)

## 保存数据
save(sce.all, file = "./1-QC/sce.all_raw.Rdata")


# 02数据质控
## 设置环境
rm(list = ls()) 
options(stringsAsFactors = F) 
library(Seurat)
library(ggplot2)
library(clustree)
library(cowplot)
library(dplyr)
library(stringr)
library(data.table)
library(ggsci)
load(file = "./1-QC/sce.all_raw.Rdata")
load(file = "./0-Rawdata/my_color.Rdata")
color1 = col_vector[1:24]

## 计算线粒体基因比例
mito_genes = rownames(sce.all)[grep("^MT-", rownames(sce.all))] # 人和鼠的基因名字稍微不一样 
mito_genes #13个线粒体基因
sce.all = PercentageFeatureSet(sce.all, "^MT-", col.name = "percent_mito")

## 计算核糖体基因比例
ribo_genes = rownames(sce.all)[grep("^RP[SL]", rownames(sce.all))]
ribo_genes
sce.all = PercentageFeatureSet(sce.all, "^RP[SL]", col.name = "percent_ribo")

## 计算红血细胞基因比例
rownames(sce.all)[grep("^HB[^(P)]", rownames(sce.all))]
sce.all = PercentageFeatureSet(sce.all, "^HB[^(P)]", col.name = "percent_hb")

## 可视化细胞的上述比例情况
feats <- c("nFeature_RNA", "nCount_RNA")
p1 = VlnPlot(sce.all, group.by = "orig.ident", features = feats, pt.size = 0, ncol = 2, cols = color1) + 
  NoLegend()
p1
ggsave(p1, filename="./1-QC/1-raw.nfeature.pdf", width = 12, height = 4)

feats <- c("percent_mito", "percent_ribo", "percent_hb")
p2 = VlnPlot(sce.all, group.by = "orig.ident", features = feats, pt.size = 0, ncol = 3, same.y.lims = T, cols = color1) + 
  # scale_y_continuous(breaks=seq(0, 100, 5)) +
  NoLegend()
p2	
ggsave(p2, filename="./1-QC/2-raw.pecentage.pdf", width = 18, height = 4)

p3 = FeatureScatter(sce.all, "nCount_RNA", "nFeature_RNA", group.by = "orig.ident", pt.size = 1, cols = color1)
p3
ggsave(p3, filename="./1-QC/3-raw.Scatterplot.pdf", width = 5.5, height = 3.75)


# 根据上述指标，过滤低质量细胞/基因
## 过滤指标1:最少表达基因数的细胞&最少表达细胞数的基因
selected_c <- WhichCells(sce.all, expression = nFeature_RNA > 500)
selected_f <- rownames(sce.all)[Matrix::rowSums(sce.all@assays$RNA@counts > 0) > 3]

sce.all.filt <- subset(sce.all, features = selected_f, cells = selected_c)
dim(sce.all) 
dim(sce.all.filt) 

table(sce.all@meta.data$orig.ident) 
table(sce.all.filt@meta.data$orig.ident) 
## 可以看到，主要是过滤了基因，其次才是细胞

## 过滤指标2:线粒体/核糖体基因比例(根据上面的violin图)
selected_mito <- WhichCells(sce.all.filt, expression = percent_mito < 15)
selected_ribo <- WhichCells(sce.all.filt, expression = percent_ribo > 3)
selected_hb <- WhichCells(sce.all.filt, expression = percent_hb < 1)
length(selected_hb)
length(selected_ribo)
length(selected_mito)

sce.all.filt <- subset(sce.all.filt, cells = selected_mito)
sce.all.filt <- subset(sce.all.filt, cells = selected_ribo)
sce.all.filt <- subset(sce.all.filt, cells = selected_hb)
dim(sce.all.filt)

table(sce.all.filt$orig.ident) 

## 过滤指标3:过滤特定基因
## 过滤线粒体基因
dim(sce.all.filt)
sce.all.filt <- sce.all.filt[!grepl("^MT-", rownames(sce.all.filt), ignore.case = T), ]
sce.all.filt <- sce.all.filt[!grepl("^RP[SL]", rownames(sce.all.filt), ignore.case = T), ]
dim(sce.all.filt) 
## 当然，还可以过滤更多


## 可视化过滤后的情况
feats <- c("nFeature_RNA", "nCount_RNA")
p1_filtered = VlnPlot(sce.all.filt, group.by = "orig.ident", features = feats, pt.size = 0, ncol = 2, cols = color1) + 
  NoLegend()
p1_filtered
ggsave(p1_filtered, filename="./1-QC/4-filter.nfeature.pdf", width = 12, height = 4)

feats <- c("percent_mito", "percent_ribo", "percent_hb")
p2_filtered = VlnPlot(sce.all.filt, group.by = "orig.ident", features = feats, pt.size = 0, ncol = 3, cols = color1) + 
  NoLegend()
p2_filtered
ggsave(p2_filtered, filename="./1-QC/5-filter.percentage.pdf", width = 18, height = 4)


## 细胞周期评分
sce.all.filt = NormalizeData(sce.all.filt)
s.genes = Seurat::cc.genes.updated.2019$s.genes
g2m.genes = Seurat::cc.genes.updated.2019$g2m.genes
sce.all.filt = CellCycleScoring(object = sce.all.filt, 
                                s.features = s.genes, 
                                g2m.features = g2m.genes, 
                                set.ident = TRUE)
p4 = VlnPlot(sce.all.filt, features = c("S.Score", "G2M.Score"), group.by = "orig.ident", 
             ncol = 2, pt.size = 0, cols = color1)
p4
ggsave(p4, filename = "./1-QC/6-filter.cellcycle.pdf", width = 12, height = 4)

sce.all.filt@meta.data  %>% ggplot(aes(S.Score,G2M.Score)) + geom_point(aes(color=Phase)) + theme_bw() + scale_color_nejm()
ggsave(filename = "./1-QC/7-filter.cycle.details.pdf", width = 4, height = 3)
# S.Score较高的为S期，G2M.Score较高的为G2M期，都比较低的为G1期

## 保存数据
save(sce.all.filt, file = "./1-QC/sce.all.filter.Rdata")



