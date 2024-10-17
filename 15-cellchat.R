
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
library(patchwork)
library(ggsci)
load(file = "./14-cellchat/mvi.cellchat.Rdata")
load(file = "./14-cellchat/nomvi.cellchat.Rdata")

# 合并
object.list = list(nomvi = nomvi.cellchat, mvi = mvi.cellchat)
cellchat = mergeCellChat(object.list, add.names = names(object.list))
cellchat

# 比较细胞通讯相互作用的总数和相互作用的强度
gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1,2))
gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "weight")
gg1 + gg2
ggsave(filename = "./14-cellchat/1-compare.pdf", width = 3.5, height = 2.6)

# 不同细胞群间相互作用次数或相互作用强度的差异
## 以circle plot的形式展示第二个组别中相较于第一个组别细胞通讯发生的变化，红色为上调蓝色为下调
pdf("./14-cellchat/2-net_number_strength.pdf")
par(mfrow = c(1,2), xpd=TRUE)
netVisual_diffInteraction(cellchat, weight.scale = T)
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight")
dev.off()

## 热图更详细地显示相互作用的差异数量或相互作用强度
gg1 <- netVisual_heatmap(cellchat,comparison = c(1,2))
gg2 <- netVisual_heatmap(cellchat, measure = "weight",comparison = c(1,2))
gg1 + gg2
ggsave(filename = "./14-cellchat/3-compare-heatmap.pdf", width = 5.9, height = 3)

## 分开展示各组的通讯情况
weight.max <- getMaxWeight(object.list, attribute = c("idents","count"))
pdf("./14-cellchat/4-net_number_strength-sigle.pdf")
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_circle(object.list[[i]]@net$count, weight.scale = T, label.edge= F, edge.weight.max = weight.max[2], edge.width.max = 12, title.name = paste0("Number of interactions - ", names(object.list)[i]))
}
dev.off()

# 不同细胞间相互作用次数或相互作用强度的差异（以Fib和B细胞为例）
group.cellType = c("Fibroblast", "Macrophage", "B cell", "T cell")
group.cellType = factor(group.cellType, levels = c("Fibroblast", "Macrophage", "B cell", "T cell"))
object.list2 = lapply(object.list, function(x) {mergeInteractions(x, group.cellType)})
cellchat = mergeCellChat(object.list2, add.names = names(object.list2))

## 各组别的通讯情况
weight.max <- getMaxWeight(object.list2, slot.name = c("idents", "net", "net"), attribute = c("idents","count", "count.merged"))
pdf("./14-cellchat/5-net_number_strength-immue.pdf")
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list2)) {
  netVisual_circle(object.list2[[i]]@net$count.merged,
                   weight.scale = T, label.edge= T,
                   edge.weight.max = weight.max[3],
                   edge.width.max = 12,
                   title.name = paste0("Number of interactions - ", names(object.list)[i]))
}
dev.off()

## 差异性circle plot通讯情况
pdf("./14-cellchat/6-net_number_compare-immue.pdf")
par(mfrow = c(1,2), xpd=TRUE)
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "count.merged",comparison = c(1,2))
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight.merged",comparison = c(1,2))
dev.off()

# 比较不同组别中的主要sources和targets
num.link = sapply(object.list, function(x) {rowSums(x@net$count) + colSums(x@net$count)-diag(x@net$count)})
weight.MinMax = c(min(num.link), max(num.link)) # control the dot size in the different datasets
gg = list()
for (i in 1:length(object.list)) {
  object.list[[i]] = netAnalysis_computeCentrality(object.list[[i]])
  gg[[i]] = netAnalysis_signalingRole_scatter(object.list[[i]], title = names(object.list)[i], weight.MinMax = weight.MinMax)
}
patchwork::wrap_plots(plots = gg)

## 具体信号变化
gg1 = netAnalysis_signalingChanges_scatter(object.list, idents.use = "Macrophage")
gg2 = netAnalysis_signalingChanges_scatter(object.list, idents.use = "B cell")
gg3 = netAnalysis_signalingChanges_scatter(object.list, idents.use = "T cell")
patchwork::wrap_plots(plots = list(gg1, gg2, gg3))


# 比较各个信号通路的总体信息流
gg1 = rankNet(cellchat, mode = "comparison", stacked = T, do.stat = TRUE,comparison = c(1,2))
gg2 = rankNet(cellchat, mode = "comparison", stacked = F, do.stat = TRUE,comparison = c(1,2))
gg1 + gg2
ggsave(filename = "./14-cellchat/9-compare-details.pdf", width = 7, height = 4.7)

# 比较与每个细胞群相关的传出outgoing（或传入incoming）信号
library(ComplexHeatmap)
i = 1
pathway.union <- union(object.list[[i]]@netP$pathways, object.list[[i+1]]@netP$pathways)

## outgoing
ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i], width = 5, height = 6)
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i+1], width = 5, height = 6)
draw(ht1 + ht2)

## incoming
ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "incoming", 
                                        signaling = pathway.union, title = names(object.list)[i],
                                        width = 5, height = 6, color.heatmap = "GnBu")
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "incoming",
                                        signaling = pathway.union, title = names(object.list)[i+1],
                                        width = 5, height = 6, color.heatmap = "GnBu")
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))

## all
ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "all", 
                                        signaling = pathway.union, title = names(object.list)[i], 
                                        width = 5, height = 6, color.heatmap = "OrRd")
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "all",
                                        signaling = pathway.union, title = names(object.list)[i+1], 
                                        width = 5, height = 6, color.heatmap = "OrRd")
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))


# 识别上调和下调的信号配体受体对
netVisual_bubble(cellchat, sources.use = 5, targets.use = c(7:9),  comparison = c(1, 2), angle.x = 45)

## 分别展示上调下调通路
gg1 = netVisual_bubble(cellchat, sources.use = 5, targets.use = c(1,2,4),  comparison = c(1, 2), max.dataset = 2, title.name = "Increased signaling in mvi", angle.x = 45, remove.isolate = T)
gg2 = netVisual_bubble(cellchat, sources.use = 5, targets.use = c(1,2,4),  comparison = c(1, 2), max.dataset = 1, title.name = "Decreased signaling in mvi", angle.x = 45, remove.isolate = T)
gg1 + gg2

# 通过差异表达分析识别功能失调信号

# define a positive dataset, i.e., the dataset with positive fold change against the other dataset
pos.dataset = "mvi"

# define a char name used for storing the results of differential expression analysis
features.name = pos.dataset

# perform differential expression analysis
cellchat <- identifyOverExpressedGenes(cellchat, 
                                       group.dataset = "datasets", 
                                       pos.dataset = pos.dataset,
                                       features.name = features.name,
                                       only.pos = FALSE,
                                       thresh.pc = 0.1,
                                       thresh.fc = 0.1,
                                       thresh.p = 1)

# map the results of differential expression analysis onto the inferred cell-cell communications to easily manage/subset the ligand-receptor pairs of interest
net <- netMappingDEG(cellchat, features.name = features.name)

# extract the ligand-receptor pairs with upregulated ligands in mvi
net.up <- subsetCommunication(cellchat, net = net,
                              datasets = "mvi",ligand.logFC = 0.2,
                              receptor.logFC = NULL)

# extract the ligand-receptor pairs with upregulated ligands and upregulated recetptors in nomvi, i.e.,downregulated in mvi
net.down <- subsetCommunication(cellchat, net = net,
                                datasets = "nomvi",
                                ligand.logFC = -0.1,
                                receptor.logFC = -0.1)

gene.up <- extractGeneSubsetFromPair(net.up, cellchat)
gene.down <- extractGeneSubsetFromPair(net.down, cellchat)


pairLR.use.up = net.up[,"interaction_name", drop = F]
gg1 <- netVisual_bubble(cellchat,
                        pairLR.use = pairLR.use.up,
                        sources.use = 5, targets.use = c(7:9),
                        comparison = c(1, 2),
                        angle.x = 90,
                        remove.isolate = T,
                        title.name = paste0("Up-regulated signaling in ", names(object.list)[2]))

pairLR.use.down = net.down[, "interaction_name", drop = F]
gg2 <- netVisual_bubble(cellchat,
                        pairLR.use = pairLR.use.down,
                        sources.use = 5, targets.use = c(7:9),
                        comparison = c(1, 2),
                        angle.x = 90, remove.isolate = T,
                        title.name = paste0("Down-regulated signaling in ", names(object.list)[2]))
gg1 + gg2

# 使用和弦图可视化上调和下调的信号配体受体对:
## Chord diagram
par(mfrow = c(1,2), xpd=TRUE)
netVisual_chord_gene(object.list[[2]],
                     sources.use = 5,targets.use = c(7:9),
                     slot.name = 'net',net = net.up, 
                     lab.cex = 0.8,small.gap = 3.5,
                     title.name = paste0("Up-regulated signaling in ",names(object.list)[2]))
netVisual_chord_gene(object.list[[1]],
                     sources.use = 5, targets.use = c(7:9),
                     slot.name = 'net', net = net.down,
                     lab.cex = 0.8, small.gap = 3.5,
                     title.name = paste0("Down-regulated signaling in ", names(object.list)[2]))

# 使用层次图、圆圈图或和弦图直观地比较细胞间的通讯

## heatmap1
pathways.show = c("ANNEXIN") 
par(mfrow = c(1,2), xpd=TRUE)
ht <- list()
for (i in 1:length(object.list)) {
  ht[[i]] <- netVisual_heatmap(object.list[[i]], signaling = pathways.show, color.heatmap = "Reds",title.name = paste(pathways.show, "signaling ",names(object.list)[i]))
}
ComplexHeatmap::draw(ht[[1]] + ht[[2]], ht_gap = unit(0.5, "cm"))

## heatmap2
pathways.show = c("CXCL") 
par(mfrow = c(1,2), xpd=TRUE)
ht <- list()
for (i in 1:length(object.list)) {
  ht[[i]] <- netVisual_heatmap(object.list[[i]], signaling = pathways.show, color.heatmap = "Reds",title.name = paste(pathways.show, "signaling ",names(object.list)[i]))
}
ComplexHeatmap::draw(ht[[1]] + ht[[2]], ht_gap = unit(0.5, "cm"))
