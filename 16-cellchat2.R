

# 设置环境
rm(list=ls())
options(stringsAsFactors = F)
library(Seurat)
library(ggplot2)
library(clustree)
library(cowplot)
library(dplyr)
library(ggpubr)
library(CellChat)

# fib数据整理
load(file = "./8-fibsub/fibsub.Rdata")
fibsub@meta.data$celltype = paste0("Fibroblast", fibsub@meta.data$seurat_clusters)

# epi数据整理
load(file = "./6-episub/episub.Rdata")
episub@meta.data$celltype = "Ductal"

# 数据合并
cell.use = merge(x = fibsub, y = episub)
rm(list = c("fibsub", "episub"))
cell.use = SetIdent(cell.use, value = "celltype")

# 输入数据
data.input = cell.use@assays$RNA@data
identity = data.frame(group = cell.use$celltype, row.names = names(cell.use$celltype)) # create a dataframe consisting of the cell labels
unique(identity$group) # check the cell labels

# 创建cellchat对象
cellchat = createCellChat(object = data.input)
cellchat

# 加入metadata
cellchat = addMeta(cellchat, meta = identity, meta.name = "labels")
cellchat = setIdent(cellchat, ident.use = "labels") # set "labels" as default cell identity
levels(cellchat@idents) # show factor levels of the cell labels
groupSize = as.numeric(table(cellchat@idents)) # number of cells in each cell group

# 导入配体受体数据库
CellChatDB = CellChatDB.human
showDatabaseCategory(CellChatDB)
colnames(CellChatDB$interaction)
dplyr::glimpse(CellChatDB$interaction)

CellChatDB.use = subsetDB(CellChatDB, search = "Secreted Signaling") # use Secreted Signaling for cell-cell communication analysis
cellchat@DB = CellChatDB.use # set the used database in the object
unique(CellChatDB$interaction$annotation)

# 预处理
cellchat = subsetData(cellchat) # subset the expression data of signaling genes for saving computation cost
future::plan("multiprocess", workers = 10) 
cellchat = identifyOverExpressedGenes(cellchat) #寻找高表达的基因#
cellchat = identifyOverExpressedInteractions(cellchat) #寻找高表达的通路
cellchat = projectData(cellchat, PPI.human) #投影到PPI


# 相互作用推断
cellchat = computeCommunProb(cellchat) #默认cutoff的值为20%，即表达比例在25%以下的基因会被认为是0， trim = 0.1可以调整比例阈值
cellchat = computeCommunProbPathway(cellchat)
cellchat = filterCommunication(cellchat, min.cells = 10)
df.netp = subsetCommunication(cellchat, slot.name = "netP")
cellchat = aggregateNet(cellchat)

# 保存数据
fibepi.cellchat = cellchat
save(fibepi.cellchat, file = "./15-cellchat/fibepi.cellchat.Rdata")


# 设置环境
rm(list=ls())
options(stringsAsFactors = F)
library(Seurat)
library(ggplot2)
library(clustree)
library(cowplot)
library(dplyr)
library(ggpubr)
library(CellChat)
load(file = "./15-cellchat/fibepi.cellchat.Rdata")
fibepi.cellchat@netP$pathways
head(fibepi.cellchat@LR$LRsig)


# 可视化
## 总体情况热图
fibepi.cellchat = netAnalysis_computeCentrality(fibepi.cellchat, slot.name = "netP")
ht1 = netAnalysis_signalingRole_heatmap(fibepi.cellchat, pattern = "outgoing")
ht2 = netAnalysis_signalingRole_heatmap(fibepi.cellchat, pattern = "incoming")
ht1 + ht2

## 总体情况气泡图
netVisual_bubble(fibepi.cellchat, sources.use = 2:3, targets.use = 1, remove.isolate = FALSE)

## 总体情况点图
gg1 = netAnalysis_signalingRole_scatter(fibepi.cellchat)
gg2 = netAnalysis_signalingRole_scatter(fibepi.cellchat, signaling = c("ANGPTL", "MK", "SPP1", "PERIOSTIN"))
gg1 + gg2

## 细胞通讯图
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(fibepi.cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(fibepi.cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")

## 具体信号通路

## Heatmap
pathways.show = c("MK") 
par(mfrow=c(1,5))
gg1 = netVisual_heatmap(fibepi.cellchat, signaling = pathways.show, color.heatmap = "Reds")

pathways.show = c("SPP1") 
par(mfrow=c(1,5))
gg2 = netVisual_heatmap(fibepi.cellchat, signaling = pathways.show, color.heatmap = "Reds")

pathways.show = c("MIF") 
par(mfrow=c(1,5))
gg3 = netVisual_heatmap(fibepi.cellchat, signaling = pathways.show, color.heatmap = "Reds")

pathways.show = c("ANGPTL") 
par(mfrow=c(1,5))
gg4 = netVisual_heatmap(fibepi.cellchat, signaling = pathways.show, color.heatmap = "Reds")

pathways.show = c("PERIOSTIN") 
par(mfrow=c(1,5))
gg5 = netVisual_heatmap(fibepi.cellchat, signaling = pathways.show, color.heatmap = "Reds")

gg1 + gg2 + gg3 + gg4 + gg5








