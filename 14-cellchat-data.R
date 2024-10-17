
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
load("./4-celltype/sce.har.db.Rdata")

# 取非转移子集
nomvi.sce = subset(x = sce.har.db, MVI == "None_MVI")
rm(list = c("sce.har.db"))

# 导入数据
nomvi.sce = SetIdent(nomvi.sce, value = "celltype")
DimPlot(nomvi.sce, reduction = "umap", group.by = "celltype", label = T)

# 输入数据
data.input = nomvi.sce@assays$RNA@data
identity = data.frame(group = nomvi.sce$celltype, row.names = names(nomvi.sce$celltype)) # create a dataframe consisting of the cell labels
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
nomvi.cellchat = cellchat
save(nomvi.cellchat, file = "./14-cellchat/nomvi.cellchat.Rdata")

# 设置环境
rm(list=ls())
load("./4-celltype/sce.har.db.Rdata")

# 取转移子集
mvi.sce = subset(x = sce.har.db, MVI == "MVI")
rm(list = c("sce.har.db"))

# 导入数据
mvi.sce = SetIdent(mvi.sce, value = "celltype")
DimPlot(mvi.sce, reduction = "umap", group.by = "celltype", label = T)

# 输入数据
data.input = mvi.sce@assays$RNA@data
identity = data.frame(group = mvi.sce$celltype, row.names = names(mvi.sce$celltype)) # create a dataframe consisting of the cell labels
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
mvi.cellchat = cellchat
save(mvi.cellchat, file = "./14-cellchat/mvi.cellchat.Rdata")
