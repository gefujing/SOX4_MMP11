
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
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggrepel)
library(enrichplot)
library(SCENIC)
library(doParallel)

load("./8-fibsub/fibsub.Rdata")
load(file = "./0-Rawdata/my_color.Rdata")
sce.markers = read.csv(file = "./8-fibsub/fibsub.markers.csv", header = T)

# 设置核心
library(future) # check the current active plan
plan()
plan("multiprocess", workers = 10)
plan()


# 准备数据库
# 保证 cisTarget_databases 文件夹下面有下载好2个1G的文件
scenicOptions <- initializeScenic(org = "hgnc", dbDir = "./10-fibsub-scenic/cisTarget_databases", nCores=10)
saveRDS(scenicOptions, file = "./10-fibsub-scenic/scenicOptions.Rds") 

# 准备输入数据
exprMat = as.matrix(fibsub@assays$RNA@data)
cellInfo = fibsub@meta.data[ ,c(1,4,5,6,8)]

### Co-expression network-耗时间
genesKept <- geneFiltering(exprMat, scenicOptions)
exprMat_filtered <- exprMat[genesKept, ]
exprMat_filtered[1:4,1:4]
dim(exprMat_filtered)
runCorrelation(exprMat_filtered, scenicOptions)
exprMat_filtered_log <- log2(exprMat_filtered + 1) 
runGenie3(exprMat_filtered_log, scenicOptions)

### Build and score the GRN
exprMat_log <- log2(exprMat+1)
scenicOptions@settings$dbs <- scenicOptions@settings$dbs["10kb"] # Toy run settings
scenicOptions <- runSCENIC_1_coexNetwork2modules(scenicOptions)
scenicOptions <- runSCENIC_2_createRegulons(scenicOptions, coexMethod = c("top5perTarget")) # Toy run settings
scenicOptions <- runSCENIC_3_scoreCells(scenicOptions, exprMat_log) 
scenicOptions <- runSCENIC_4_aucell_binarize(scenicOptions)
tsneAUC(scenicOptions, aucType = "AUC") # choose settings


# 结果解读
rm(list = ls()) 
library(Seurat) 
library(SCENIC)
library(doParallel)

scenicOptions = readRDS(file = "./int/scenicOptions.Rds")
motifEnrichment_selfMotifs_wGenes <- loadInt(scenicOptions, "motifEnrichment_selfMotifs_wGenes") 
as.data.frame(sort(table(motifEnrichment_selfMotifs_wGenes$highlightedTFs), decreasing = T)) 

# 可视化IRF7基因的motif序列特征
tableSubset <- motifEnrichment_selfMotifs_wGenes[highlightedTFs == "IRF7"]
viewMotifs(tableSubset) 

regulonTargetsInfo <- loadInt(scenicOptions, "regulonTargetsInfo")
tableSubset <- regulonTargetsInfo[TF=="IRF7" & highConfAnnot==TRUE]
viewMotifs(tableSubset) 

# 作图数据
rm(list = ls())
library(Seurat) 
library(SCENIC)

load("./8-fibsub/fibsub.Rdata")
load(file = "./0-Rawdata/my_color.Rdata")
tfmat = readRDS(file = "./10-fibsub-scenic/int/3.4_regulonAUC.Rds")
tfmat = as.data.frame(tfmat@assays@data@listData[["AUC"]])
identical(colnames(tfmat), colnames(fibsub))
colanno = fibsub@meta.data[,c(6,8)]
k = str_detect(string = rownames(tfmat), pattern = "extended")
tfmat = tfmat[!k,]

## 差异分析

library(limma)
Group = colanno$seurat_clusters
Group = factor(Group, levels = c(0,1))
design = model.matrix(~Group)
fit = lmFit(tfmat, design)
fit = eBayes(fit)
deg = topTable(fit, coef=2, number = Inf)

k = c("NUAK1 (89g)", "SOX4 (36g)", "MAF (25g)", "CREB3L1 (55g)", "UQCRB (12g)",
      "CEBPD (27g)", "FOSB (129g)", "ATF3 (10g)", "CREM (72g)", "STAT1 (80g)")

## 火山图
tfmat = tfmat[k,]
rownames(tfmat) = c("NUAK1", "SOX4", "MAF", "CREB3L1", "UQCRB", "CEBPD", "FOSB", "ATF3", "CREM", "STAT1")
colanno = colanno[order(colanno$seurat_clusters),]
tfmat = tfmat[,rownames(colanno)]

dat.exp = t(scale(t(tfmat)))
dat.exp[dat.exp>=2] = 2
dat.exp[dat.exp<=(-2)] = -2

pheatmap(dat.exp,
         color = colorRampPalette(c("#42218E", "#7B3B8C", "#F9CE00"))(99),
         cluster_rows = T,
         cluster_cols = F,
         show_colnames = F,
         annotation_col = colanno
)


# dimplot
rm(list = ls())
library(Seurat) 
library(SCENIC)
library(patchwork)

load("./8-fibsub/fibsub.Rdata")
load(file = "./0-Rawdata/my_color.Rdata")
tfmat = readRDS(file = "./10-fibsub-scenic/int/3.4_regulonAUC.Rds")
tfmat = as.data.frame(tfmat@assays@data@listData[["AUC"]])

tsne = as.data.frame(fibsub@reductions[["tsne"]]@cell.embeddings)
identical(colnames(tfmat), rownames(tsne))
k = c("NUAK1 (89g)", "SOX4 (36g)", "MAF (25g)", "CREB3L1 (55g)", "UQCRB (12g)",
      "CEBPD (27g)", "FOSB (129g)", "ATF3 (10g)", "CREM (72g)", "STAT1 (80g)")

tfmat = tfmat[k,]
tfmat = t(tfmat)
colnames(tfmat) = c("NUAK1", "SOX4", "MAF", "CREB3L1", "UQCRB", "CEBPD", "FOSB", "ATF3", "CREM", "STAT1")

dimdat = cbind(tsne, tfmat)

p1 = ggplot(data = dimdat, mapping = aes(x = tSNE_1, y = tSNE_2, color = CREB3L1))+
  geom_point(size = 0.5)+
  theme_bw()+
  scale_color_material("red")

p2 = ggplot(data = dimdat, mapping = aes(x = tSNE_1, y = tSNE_2, color = MAF))+
  geom_point(size = 0.5)+
  theme_bw()+
  scale_color_material("red")

p3 = ggplot(data = dimdat, mapping = aes(x = tSNE_1, y = tSNE_2, color = NUAK1))+
  geom_point(size = 0.5)+
  theme_bw()+
  scale_color_material("red")

p4 = ggplot(data = dimdat, mapping = aes(x = tSNE_1, y = tSNE_2, color = SOX4))+
  geom_point(size = 0.5)+
  theme_bw()+
  scale_color_material("red")

p5 = ggplot(data = dimdat, mapping = aes(x = tSNE_1, y = tSNE_2, color = UQCRB))+
  geom_point(size = 0.5)+
  theme_bw()+
  scale_color_material("red")

pall = p1 + p2 +p3 +p4 +p5
ggsave(pall, filename = "./10-fibsub-scenic/3-c1-tfs.pdf", width = 8.86, height = 4.16)

pall = p3 + p4
ggsave(pall, filename = "./10-fibsub-scenic/4-c1-tfs-2.pdf", width = 6.83, height = 2.58)


p1 = ggplot(data = dimdat, mapping = aes(x = tSNE_1, y = tSNE_2, color = CEBPD))+
  geom_point(size = 0.5)+
  theme_bw()+
  scale_color_material("blue")

p2 = ggplot(data = dimdat, mapping = aes(x = tSNE_1, y = tSNE_2, color = STAT1))+
  geom_point(size = 0.5)+
  theme_bw()+
  scale_color_material("blue")

p3 = ggplot(data = dimdat, mapping = aes(x = tSNE_1, y = tSNE_2, color = ATF3))+
  geom_point(size = 0.5)+
  theme_bw()+
  scale_color_material("blue")

p4 = ggplot(data = dimdat, mapping = aes(x = tSNE_1, y = tSNE_2, color = FOSB))+
  geom_point(size = 0.5)+
  theme_bw()+
  scale_color_material("blue")

p5 = ggplot(data = dimdat, mapping = aes(x = tSNE_1, y = tSNE_2, color = CREM))+
  geom_point(size = 0.5)+
  theme_bw()+
  scale_color_material("blue")

pall = p1 + p2 +p3 +p4 +p5
ggsave(pall, filename = "./10-fibsub-scenic/5-c0-tfs.pdf", width = 8.86, height = 4.16)

pall = p1 + p3
ggsave(pall, filename = "./10-fibsub-scenic/6-c0-tfs-2.pdf", width = 6.83, height = 2.58)
  



