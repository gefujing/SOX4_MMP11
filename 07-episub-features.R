
# 准备环境
rm(list=ls()) 
options(stringsAsFactors = F) 
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
library(GSVA)
library(msigdbr)
load("./6-episub/episub.Rdata")
load(file = "./0-Rawdata/my_color.Rdata")
sce.markers = read.csv(file = "./6-episub/episub.markers.csv", header = T)

# 功能聚类(基于marker gene)
## clusterprofiler
s2e = bitr(sce.markers$gene,'SYMBOL','ENTREZID','org.Hs.eg.db') ## 将SYMBOL转成ENTREZID
sce.markers = merge(sce.markers, s2e, by.x='gene', by.y='SYMBOL')

## 函数split()可以按照分组因子，把向量，矩阵和数据框进行适当的分组。
## 它的返回值是一个列表，代表分组变量每个水平的观测。
gcSample = split(sce.markers$ENTREZID, sce.markers$cluster) 

## KEGG
kk <- compareCluster(gcSample, fun = "enrichKEGG", organism = "hsa", pvalueCutoff = 0.05)
p <- dotplot(kk)
ggsave(p, filename = "./6-episub/6-cluster-kegg.pdf", width = 5.9, height = 7.8)

kk <- compareCluster(gcSample, fun = "enrichGO", OrgDb = "org.Hs.eg.db", ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 0.01, qvalueCutoff = 0.05)
p <- dotplot(kk)
ggsave(p, filename = "./6-episub/6-cluster-go.pdf", width = 5.9, height = 7.5)


## 可视化
temp = kk@compareClusterResult
fz = as.numeric(str_split(string = temp$GeneRatio, pattern = "/", simplify = T)[,1])
fm = as.numeric(str_split(string = temp$GeneRatio, pattern = "/", simplify = T)[,2])
gr = fz/fm
fz = as.numeric(str_split(string = temp$BgRatio, pattern = "/", simplify = T)[,1])
fm = as.numeric(str_split(string = temp$BgRatio, pattern = "/", simplify = T)[,2])
br = fz/fm
ge = gr/br
temp$generatio = gr
temp$enrichment = ge

top3 = temp %>% group_by(Cluster) %>% top_n(n=3, wt = p.adjust)
top3$Description = factor(top3$Description , levels = unique(top3$Description))
top3 = top3[c(-24, -19, -23),]

ggplot(data = top3, mapping = aes(x = Cluster, y = Description, size = generatio, color = p.adjust))+
  geom_point()+
  scale_color_gradientn(colours = rev(c("#FFD92F","#FEE391",brewer.pal(11, "Spectral")[7:11]))) +
  theme_bw()
ggsave(filename = "./6-episub/6-cluster-go2.pdf", width = 6.5, height = 4.7)  


# 功能聚类（基于deg）
## 找差异基因
deg = FindMarkers(object = episub, ident.1 = c(2,4), ident.2 = c(0,1,3,5,6), test.use='MAST' )  ## MAST在单细胞领域较为常用
VlnPlot(episub, features = c("CDH1", "OCLN", "VIM", "CTNNB1", "S100A4", "FN1"), pt.size =  0)
FeaturePlot(episub, reduction = "tsne", features = c("FN1"), cols = c("grey", "red"))
ggsave(filename = "./6-episub/6-tsne-emt.pdf", width = 3, height = 2.7) 
VlnPlot(episub, features = c("FN1"), pt.size = 0) + scale_fill_manual(values = col_vector)
ggsave(filename = "./6-episub/6-vln-emt.pdf", width = 3.3, height = 2.7) 
save(deg, file = "./6-episub/DEGs-for-24-01356-cluster.Rdata")

## volcano
degdf = deg
degdf$symbol = rownames(deg)
logFC_t = 1
P.Value_t = 1e-100
degdf$change = ifelse(degdf$p_val_adj < P.Value_t & degdf$avg_log2FC < -logFC_t, "Down",
                      ifelse(degdf$p_val_adj < P.Value_t & degdf$avg_log2FC > logFC_t,"Up","Stable"))
degdf$perchange = deg$pct.1 - deg$pct.2

degdf$label = ""
degdf = degdf[order(abs(degdf$avg_log2FC), decreasing = T), ]
up.genes = head(degdf$symbol[which(degdf$change == "Up")],5)
down.genes = head(degdf$symbol[which(degdf$change == "Down")],5)
top10genes = c(as.character(up.genes), as.character(down.genes))
degdf$label[match(top10genes,degdf$symbol)] <- top10genes

ggplot(data = degdf, mapping = aes(x = avg_log2FC, y = perchange, color = change))+
  geom_point()+
  geom_text_repel(aes(label = label), max.overlaps = Inf, alpha = 0.8)+
  theme_bw()+
  scale_color_manual(values = c("darkblue", "grey", "darkred"))

ggsave(filename = "./6-episub/7-DEGs-for-24-01356-cluster.pdf", width = 4.5, height = 3)

## KEGG
## 获取上下调基因
gene_up = rownames(deg[deg$avg_log2FC > 0,])
gene_down = rownames(deg[deg$avg_log2FC < 0,])
## 把SYMBOL改为ENTREZID
gene_up = as.character(na.omit(AnnotationDbi::select(org.Hs.eg.db, keys = gene_up, columns = 'ENTREZID', keytype = 'SYMBOL')[,2]))
gene_down = as.character(na.omit(AnnotationDbi::select(org.Hs.eg.db, keys = gene_down, columns = 'ENTREZID', keytype = 'SYMBOL')[,2]))
gene_diff = c(gene_up, gene_down)
## KEGG
gene_diff = unique(gene_diff)
kk.diff = enrichKEGG(gene = gene_diff, organism = "hsa", pvalueCutoff = 0.9, qvalueCutoff = 0.9)
dotplot(kk.diff)
ggsave(filename = "./6-episub/8-DEGs-KEGG.pdf", width = 5.3, height = 4.3)

## GO
go.diff <- enrichGO(gene = gene_diff, OrgDb = org.Hs.eg.db, ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 0.99, qvalueCutoff = 0.99, readabl = TRUE)
barplot(go.diff)
ggsave(filename = "./6-episub/9-DEGs-GO.pdf", width = 5, height = 4.5)

##GSEA

## 上一步差异分析得到差异基因列表deg后取出，p值和log2FC
nrDEG = deg[ ,c("avg_log2FC", "p_val")]
colnames(nrDEG) = c("log2FoldChange", "pvalue") ##更改列名

## 把SYMBOL转换为ENTREZID，可能有部分丢失
gene <- bitr(rownames(nrDEG), fromType = "SYMBOL", toType =  "ENTREZID", OrgDb = org.Hs.eg.db)

## 基因名、ENTREZID、logFC一一对应起来
gene$logFC = nrDEG$log2FoldChange[match(gene$SYMBOL, rownames(nrDEG))]

## 构建genelist
geneList = gene$logFC
names(geneList) = gene$ENTREZID 
geneList = sort(geneList, decreasing = T) # 降序，按照logFC的值来排序

## GSEA分析
### kegg
kk_gse = gseKEGG(geneList = geneList, organism = 'hsa', nPerm = 1000, minGSSize = 10, pvalueCutoff = 1, verbose = FALSE)
kk_gse = DOSE::setReadable(kk_gse, OrgDb='org.Hs.eg.db', keyType='ENTREZID')
sortkk = kk_gse[order(kk_gse$enrichmentScore, decreasing = T),]
setID = c("hsa04668")

gseaplot2(x = kk_gse, geneSetID = setID, pvalue_table = TRUE)
ggsave(filename = "./6-episub/10-DEGs-Gsea.pdf", width = 4.5, height = 4)

### go
kk_gse = gseGO(geneList = geneList, ont = "all", OrgDb = org.Hs.eg.db)
kk_gse = DOSE::setReadable(kk_gse, OrgDb = 'org.Hs.eg.db', keyType = 'ENTREZID')
sortkk = kk_gse[order(kk_gse$enrichmentScore, decreasing = T),]
setID = c("GO:0045087")

gseaplot2(x = kk_gse, geneSetID = setID, pvalue_table = TRUE)
ggsave(filename = "./6-episub/11-DEGs-Gsea.pdf", width = 4.5, height = 4)

## gsva
rm(list=ls()) 
library(GSVA)
library(limma)
options(stringsAsFactors = F) 
load("./6-episub/episub.Rdata")
load(file = "./0-Rawdata/my_color.Rdata")

### 准备基因集
human_KEGG = msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:KEGG") %>% dplyr::select(gs_name,gene_symbol)
human_KEGG  = human_KEGG %>% split(x = .$gene_symbol, f = .$gs_name) #后续gsva要求是list，所以将他转化为list

## 准备表达矩阵和临床信息
exp = AverageExpression(object = episub, group.by = "orig.ident") #提取count矩阵
exp = exp[["RNA"]]
meta = episub@meta.data #分组信息，为了后续作图

## GSVA
epi.gsva = gsva(expr = exp, gset.idx.list = human_KEGG, kcdf = "Poisson", parallel.sz = 10)

## limma差异分析
meta = meta[!duplicated(meta$orig.ident),]
group = meta$MVI

design = model.matrix(~0+factor(group))
colnames(design) = levels(factor(group))
rownames(design) = colnames(epi.gsva)
contrast.matrix = makeContrasts(contrasts = paste0("MVI", '-', "None_MVI"),  #"exp/ctrl"
                                levels = design)

fit1 = lmFit(epi.gsva, design)                 #拟合模型
fit2 = contrasts.fit(fit1, contrast.matrix) #统计检验
efit = eBayes(fit2)                         #修正

summary(decideTests(efit, lfc=1, p.value = 1)) #统计查看差异结果
tempOutput = topTable(efit, coef = paste0("MVI", '-', "None_MVI"), n = Inf)
degs = na.omit(tempOutput) 
write.csv(degs, file = "./6-episub/gsva_kegg_degs.results.csv")

## 可视化
pdata = degs
pdata$ID = str_replace_all(string = rownames(pdata), pattern = "_", replacement = " ") %>% str_to_title()
pdata$ID = str_remove_all(string = pdata$ID, pattern = "Kegg") %>% str_to_title()
pdata$ID = str_remove(string = pdata$ID, pattern = " ") %>% str_to_title()
pdata = pdata[order(x = pdata$t, decreasing = T),]
pdata$rank = 1:186

up = c("Nicotinate And Nicotinamide Metabolism", "Glycolysis Gluconeogenesis",
       "Homologous Recombination", "Mismatch Repair", 
       "Pentose Phosphate Pathway", "Focal Adhesion",
       "Dna Replication", "Ecm Receptor Interaction",
       "Pathways In Cancer", "Vegf Signaling Pathway")

normal = pdata$ID[sample(c(47:146),10)]
down = pdata$ID[177:186]
pp = c(up, normal, down)

rownames(pdata) = pdata$ID
pdata = pdata[pp, ]

pdata$group = ifelse(pdata$t > 1, "UP", ifelse(pdata$t < -1, "DOWN", "STABLE"))
pdata = pdata[order(x = pdata$t, decreasing = F),]
pdata$ID = factor(pdata$ID, levels = pdata$ID)


ggplot(data = pdata, aes(x = ID, y = t, fill = group)) +
  geom_col() +
  coord_flip() + 
  scale_fill_manual(values = c("UP" = "#36638a", "STABLE" = "#cccccc", "DOWN" = "#7bcd7b")) +
  xlab("") + 
  ylab("T score") + 
  guides(fill = "none") +
  theme_bw()
ggsave(filename = "./6-episub/11-MVI-GSVA.pdf", width = 5, height = 5.3)


## AUCell
rm(list=ls()) 
library(AUCell)
library(Seurat)
library(SeuratData)
library(tidyverse)
library(patchwork)
options(stringsAsFactors = F) 
load("./6-episub/episub.Rdata")
load(file = "./0-Rawdata/my_color.Rdata")

### 准备基因集
human_KEGG = msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:KEGG") %>% dplyr::select(gs_name,gene_symbol)
human_KEGG = human_KEGG %>% split(x = .$gene_symbol, f = .$gs_name) #后续gsva要求是list，所以将他转化为list

### AUCell评分并分配活性阈值
exprMatrix = episub@assays$RNA@data
cells_rankings = AUCell_buildRankings(exprMatrix)
cells_AUC = AUCell_calcAUC(geneSets = human_KEGG, rankings = cells_rankings, nCores = 10, aucMaxRank = nrow(cells_rankings)*0.05)
temp = as.data.frame(cells_AUC@assays@data@listData[["AUC"]])

### 可视化
pdata = as.data.frame(episub@reductions[["tsne"]]@cell.embeddings)
ID = c("KEGG_ECM_RECEPTOR_INTERACTION", "KEGG_ADHERENS_JUNCTION", "KEGG_CELL_ADHESION_MOLECULES_CAMS",
       "KEGG_COMPLEMENT_AND_COAGULATION_CASCADES", "KEGG_FOCAL_ADHESION")
temp = temp[ID,]
temp = as.data.frame(t(temp))

pdata = cbind(pdata, temp)
pdata = pivot_longer(pdata, cols = 3:7)

ggplot(data = pdata,mapping = aes(x = tSNE_1, y = tSNE_2, color = value))+
  geom_point()+
  scale_color_gradientn(values = seq(0,1,0.2), colours = c('#333366',"#6666FF",'#CC3333','#FFCC33'))+
  facet_grid(~name)+
  theme_bw()


