

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
library(monocle)

load("./8-fibsub/fibsub.Rdata")
load(file = "./0-Rawdata/my_color.Rdata")


## 构建monocle对象 
## 因为monocle包是使用CellDataSet对象，所以重新构建
sample_ann <- fibsub@meta.data
head(sample_ann)

gene_ann <- data.frame(
  gene_short_name = rownames(fibsub@assays$RNA) , 
  row.names = rownames(fibsub@assays$RNA))
head(gene_ann)

## 需要三个文件phenoData/featureData/counts
pd <- new("AnnotatedDataFrame", data = sample_ann)
fd <- new("AnnotatedDataFrame", data = gene_ann)
ct <- as(as.matrix(fibsub@assays$RNA@counts), "sparseMatrix")
ct[1:4,1:4]

cds <- newCellDataSet(cellData = ct,
                      phenoData = pd,
                      featureData = fd,
                      expressionFamily = negbinomial.size())
cds  #这个CellDataSet对象一定要认识清楚，务必花两个小时去摸索它。
save(cds, col_vector, file = "./11-monocle/fibsub.monocle.Rdata")


## monocle标准流程
rm(list=ls()) 
options(stringsAsFactors = F)
library(monocle)
library(Seurat)
load(file = "./11-monocle/fibsub.monocle.Rdata")

cds <- detectGenes(cds, min_expr = 1) 
cds <- cds[fData(cds)$num_cells_expressed > 10, ] # 数值可以自行摸索
cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds) 

# 并不是所有的基因都有作用，所以先进行挑选，合适的基因用来进行聚类。
disp_table <- dispersionTable(cds)
unsup_clustering_genes <- subset(disp_table, mean_expression >= 0.1)
cds <- setOrderingFilter(cds, unsup_clustering_genes$gene_id)
plot_ordering_genes(cds) 
ggsave(filename = "./11-monocle/1-qc.pdf", width = 5, height = 3.5)
plot_pc_variance_explained(cds, return_all = F) # norm_method='log'
cds <- reduceDimension(cds, max_components = 2, num_dim = 10, reduction_method = 'tSNE', verbose = T) # 其中 num_dim 参数选择基于上面的PCA图

cds <- clusterCells(cds, num_clusters = 2) 
plot_cell_clusters(cds, 1, 2)
table(pData(cds)$Cluster) 
colnames(pData(cds)) 
table(pData(cds)$seurat_clusters)
table(pData(cds)$Cluster, pData(cds)$seurat_clusters)

# 可以看到 monocle 给细胞重新定义了亚群，亚群数量是自己选择的
# 整体来说，monocle和seurat 各自独立流程定义的亚群的一致性还不错
# 只是跑流程而已

## 保存数据
save(cds, file = "./11-monocle/monocle.input.cds.Rdata")


# 运行monocle
rm(list = ls()) 
load(file = "./11-monocle/monocle.input.cds.Rdata")
cds 

# 接下来很重要，到底是看哪个性状的轨迹
table(pData(cds)$seurat_clusters)
pData(cds)$seurat_clusters = factor(pData(cds)$seurat_clusters, levels = c(0,1))

## 我们这里并不能使用 monocle的分群
# 还是依据前面的 seurat分群, 也就是说前面的代码仅仅是流程而已，我们没有使用那些结果哦

# 其实取决于自己真实的生物学意图
pData(cds)$Cluster = pData(cds)$seurat_clusters
table(pData(cds)$Cluster)
diff_test_res <- differentialGeneTest(cds, fullModelFormulaStr = "~seurat_clusters")

# Select genes that are significant at an FDR < 10%
sig_genes <- subset(diff_test_res, qval < 0.1)
sig_genes = sig_genes[order(sig_genes$pval),]
head(sig_genes[,c("gene_short_name", "pval", "qval")] ) 
cg = as.character(head(sig_genes$gene_short_name)) 
#  挑选差异最显著的基因可视化
plot_genes_jitter(cds[cg,], grouping = "seurat_clusters", color_by = "seurat_clusters", nrow= 3, ncol = NULL)
cg2 = as.character(tail(sig_genes$gene_short_name)) 
plot_genes_jitter(cds[cg2,], grouping = "seurat_clusters", color_by = "seurat_clusters", nrow= 3, ncol = NULL )

# 前面是找差异基因，后面是做拟时序分析

# 第一步: 挑选合适的基因. 有多个方法，例如提供已知的基因集，
# 这里选取统计学显著的差异基因列表
ordering_genes <- row.names(subset(diff_test_res, qval < 0.01))
ordering_genes
cds <- setOrderingFilter(cds, ordering_genes)
plot_ordering_genes(cds) 

# 第二步: 降维。降维的目的是为了更好的展示数据。函数里提供了很多种方法,
# 不同方法的最后展示的图都不太一样, 其中“DDRTree”是Monocle2使用的默认方法
cds <- reduceDimension(cds, max_components = 2,  method = 'DDRTree')

# 第三步: 对细胞进行排序
trace("project2MST", edit = T, where = asNamespace("monocle"))
cds$State = ifelse(cds$seurat_clusters == 0, 0, 1)
cds <- orderCells(cds, root_state = 0)

# 最后两个可视化函数，对结果进行可视化
plot_cell_trajectory(cds, color_by = "seurat_clusters", size=1, show_backbone=TRUE)
ggsave(filename = "./11-monocle/2-monocle_cell_trajectory_for_seurat.pdf", width = 4.8, height = 3.5)
length(cg)
plot_genes_in_pseudotime(cds[cg,], color_by = "seurat_clusters") 
ggsave(filename = "./11-monocle/3-monocle_plot_genes_in_pseudotime_for_seurat.pdf", width = 5.5, height = 5)

phe = as.data.frame(pData(cds))
phe$seurat_clusters = factor(phe$seurat_clusters, levels = c(0, 1))
boxplot(phe$Pseudotime, phe$seurat_clusters)

# https://davetang.org/muse/2017/10/01/getting-started-monocle/
# 前面根据差异基因，推断好了拟时序，也就是说把差异基因动态化了
# 后面就可以具体推断哪些基因随着拟时序如何的变化
# pseudotime is now a column in the phenotypic data as well as the cell state
head(pData(cds))
# 这个differentialGeneTest会比较耗费时间
my_pseudotime_de <- differentialGeneTest(cds, fullModelFormulaStr = "~sm.ns(Pseudotime)", cores = 10)
# 不知道为什么在Mac电脑无法开启并行计算了 ，不过我测试了在Windows 电脑设置cores = 4是可以的
# 如果你是Mac电脑，自己修改 cores = 1 即可 
head(my_pseudotime_de)
save(cds, my_pseudotime_de, file = "./11-monocle/output_of_monocle.Rdata")


# step4:可视化
## 对前面的结果进行精雕细琢
rm(list=ls())
options(stringsAsFactors = F)
library(Seurat)
library(gplots)
library(ggplot2)
library(monocle)
library(dplyr)
library(ggsci)
load(file = "./11-monocle/output_of_monocle.Rdata")
load(file = "./0-Rawdata/my_color.Rdata")
phe = pData(cds)

# Cluster0-genes
my_pseudotime_gene = c("C7", "MGP", "APOD", "PLA2G2A", "PTGDS")
my_pseudotime_gene

plot_genes_in_pseudotime(cds[rownames(cds) %in% my_pseudotime_gene,], color_by = "Cluster") + scale_color_manual(values = col_vector)
ggsave(filename = "./11-monocle/4-c0-genes-pseudotime.pdf", width = 5.06, height = 5.19)

plot_genes_jitter(cds[rownames(cds) %in% my_pseudotime_gene,], grouping = "Cluster", color_by = "Cluster", ncol = NULL ) + scale_color_manual(values = col_vector)
ggsave(filename = "./11-monocle/5-c0-genes-pseudotime.pdf", width = 5.6, height = 4)

# Cluster1-genes
my_pseudotime_gene = c("COL11A1", "MMP11", "GJB2", "SDC1", "C1QTNF3")
my_pseudotime_gene

plot_genes_in_pseudotime(cds[rownames(cds) %in% my_pseudotime_gene,], color_by = "Cluster") + scale_color_manual(values = col_vector)
ggsave(filename = "./11-monocle/6-c1-genes-pseudotime.pdf", width = 5.06, height = 5.19)

plot_genes_jitter(cds[rownames(cds) %in% my_pseudotime_gene,], grouping = "Cluster", color_by = "Cluster", ncol = NULL ) + scale_color_manual(values = col_vector)
ggsave(filename = "./11-monocle/7-c1-genes-pseudotime.pdf", width = 5.6, height = 4)


# cluster the top 50 genes that vary as a function of pseudotime
cg = c("C7", "MGP", "APOD", "PLA2G2A", "PTGDS", "COL11A1", "MMP11", "GJB2", "SDC1", "C1QTNF3")
gene_to_cluster = my_pseudotime_de %>% arrange(qval) %>% head(100)
gene_to_cluster = gene_to_cluster$gene_short_name
gene_to_cluster
colnames(pData(cds))
table(pData(cds)$Cluster, pData(cds)$State) 
ac = pData(cds)[c("PNI","MVI","stage","Cluster",'Pseudotime')]
head(ac)
# 这个热图绘制的并不是纯粹的细胞基因表达量矩阵，而是被 Pseudotime 好了的100列，50行的矩阵
my_pseudotime_cluster <- plot_pseudotime_heatmap(cds[gene_to_cluster,],
                                                 num_clusters = 2, 
                                                # add_annotation_col = ac,
                                                 show_rownames = TRUE,
                                                 return_heatmap = TRUE)
my_pseudotime_cluster

# 聚类
rm(list=ls())
options(stringsAsFactors = F)
library(clusterProfiler)
library(org.Hs.eg.db)
library(Seurat)
library(gplots)
library(ggplot2)
library(monocle)
library(dplyr)
library(ggsci)
library(GOplot)
library(ggradar)

load(file = "./11-monocle/output_of_monocle.Rdata")
load(file = "./0-Rawdata/my_color.Rdata")
gene_to_cluster = my_pseudotime_de %>% arrange(qval) %>% head(100)
gene_to_cluster = gene_to_cluster$gene_short_name

# ID转换
P1 = c("CCL2", "FRZB", "CST3", "IGFBP7", "SFRP4", "SERPINF1", "CCDC80", "C1S", "IGFBP4",
       "SPARCL1", "CYP1B1", "MGP", "EFEMP1", "TSHZ2", "PLAC9", "CRISPLD2", "SEPP1", "C1R", "C3")
table(P1 %in% gene_to_cluster)

P2 = c("TIMP1", "SERPINE2", "APOE", "THBS4", "IGFBP6", "IGF1", "RARRES1", "APOD", "CTSC", "TGM2", "COL14A1", "A2M")
table(P2 %in% gene_to_cluster)

P3 = c("DPT", "CFH", "CXCL12", "PI16", "CFD", "CLU", "CCL19", "SRPX", "MT1M", "ADH1B", "GPX3",
       "PDK4", "HSD11B1", "SOD3", "MFAP4", "NFIB", "HSPB6", "OGN", "GSN", "PLA2G2A", "GPC3",
       "C7", "PTGDS", "BTG2", "FBLN5", "PRSS1", "NTRK2", "ABCA8", "FXYD1", "COLEC11", "ITM2A")
table(P3 %in% gene_to_cluster)

P4 = c("TGFBI", "TMEM158", "SLC16A3", "ADM", "FTH1", "GAPDH", "MIF", "RP11-400N13.3", "ENO1", "PGK1")
table(P4 %in% gene_to_cluster)  
  
P5 = c("RGS16", "CTHRC1", "COL5A2", "COL10A1","PLAU","MMP11", "THBS2", "INHBA", "COL1A1", "IGFBP3",
       "POSTN", "EPYC", "C1QTNF3", "COL12A1", "SUGCT", "IGFL2", "COL11A1", "SDC1", "MMP14", "FN1",
       "LGALS1", "LOXL2", "GREM1", "CD55", "GJB2")  
table(P5 %in% gene_to_cluster)   

P1 = bitr(P1, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
P1 = P1$ENTREZID

P2 = bitr(P2, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
P2 = P2$ENTREZID

P3 = bitr(P3, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
P3 = P3$ENTREZID

P4 = bitr(P4, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
P4 = P4$ENTREZID

P5 = bitr(P5, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
P5 = P5$ENTREZID


# 富集
k1 = enrichGO(gene = P1, OrgDb = "org.Hs.eg.db", ont = "BP", pvalueCutoff = 0.01)
k2 = enrichGO(gene = P2, OrgDb = "org.Hs.eg.db", ont = "BP", pvalueCutoff = 0.01)
k3 = enrichGO(gene = P3, OrgDb = "org.Hs.eg.db", ont = "BP", pvalueCutoff = 0.01)
k4 = enrichGO(gene = P4, OrgDb = "org.Hs.eg.db", ont = "BP", pvalueCutoff = 0.01)
k5 = enrichGO(gene = P5, OrgDb = "org.Hs.eg.db", ont = "BP", pvalueCutoff = 0.01)

# 可视化
k1 = k1@result
k2 = k2@result
k3 = k3@result
k4 = k4@result
k5 = k5@result

pid = c("GO:2000425", "GO:0006959", "GO:0010544", "GO:0045861", "GO:0006956", "GO:0002253",
        "GO:0006096", "GO:0006757", "GO:0030199", "GO:0030198")

k1 = k1[pid, ]
k2 = k2[pid, ]
k3 = k3[pid, ]
k4 = k4[pid, ]
k5 = k5[pid, ]


pdat = data.frame(row.names = pid,
                  P1 = k1$GeneRatio,
                  P2 = k2$GeneRatio,
                  P3 = k3$GeneRatio,
                  P4 = k4$GeneRatio,
                  P5 = k5$GeneRatio)

pdat = data.frame(row.names = pid,
                  P1 = str_split(string = pdat$P1, pattern = "/", simplify = T)[,1],
                  P2 = str_split(string = pdat$P2, pattern = "/", simplify = T)[,1],
                  P3 = str_split(string = pdat$P3, pattern = "/", simplify = T)[,1],
                  P4 = str_split(string = pdat$P4, pattern = "/", simplify = T)[,1],
                  P5 = str_split(string = pdat$P5, pattern = "/", simplify = T)[,1])

pdat[is.na(pdat)] = 0

pdat = as.data.frame(t(pdat))

for(i in 1:ncol(pdat)){
  pdat[,i] = as.numeric(pdat[,i])
}

pdat$group = rownames(pdat)

pdat = pdat[,c(11,1:10)]
pdat$group = factor(pdat$group, levels = c("P1", "P2", "P3", "P4", "P5"))

ggradar(pdat,
        grid.max = max(pdat[,-1]),                 # 设置坐标轴的最大值
        grid.mid = max(pdat[,-1])/2,               # 设置坐标轴的中间值
        grid.min = 0,                            # 设置坐标轴的最小值
        grid.label.size = 4,                     # 坐标轴百分比标签大小
        axis.label.size = 5,                     # 组名标签字体大小
        background.circle.colour = "white",      # 设置背景颜色
        group.point.size = 2,                    # 点大小
        group.line.width = 1,                    # 线条粗细
        plot.legend = T,                         # 是否显示图例
        legend.position = "right",               # 图例位置"top", "right", "bottom", "left"
        legend.title = "",                       # 图例标题
        legend.text.size = 10,                   # 图例文字大小
        plot.title   = "Title",                  # 标题名称
        plot.extent.x.sf = 1.2,                  # 设置图片横向延伸空间，防止外圈文字显示不全
        plot.extent.y.sf = 1.2,                  # 设置图片纵向延伸空间，防止外圈文字显示不全
) + 
  scale_color_manual(values = col_vector)
ggsave(filename = "./11-monocle/14-P1-P5-cluster.pdf", width = 9, height = 4.5)


# 亚群和拟时间
plot_cell_trajectory(cds, color_by = "seurat_clusters", size=1, show_backbone=TRUE) + scale_color_aaas()
ggsave(filename = "./11-monocle/9-monocle_cell_trajectory_for_seurat.pdf", width = 4.8, height = 3.5)

plot_cell_trajectory(cds, color_by = "Pseudotime", size=1, show_backbone=TRUE) + scale_color_continuous(type = "viridis")
ggsave(filename = "./11-monocle/10-monocle_cell_trajectory_for_Pseudotime.pdf", width = 4.8, height = 3.5)

# 单个基因
plot_cell_trajectory(cds, markers = "C7", use_color_gradient = T)
ggsave(filename = "./11-monocle/11-monocle_cell_trajectory_for_C7.pdf", width = 4.8, height = 3.5)

plot_cell_trajectory(cds, markers = "MMP11", use_color_gradient = T)
ggsave(filename = "./11-monocle/12-monocle_cell_trajectory_for_MMP11.pdf", width = 4.8, height = 3.5)

plot_cell_trajectory(cds, markers = "COL11A1", use_color_gradient = T)
ggsave(filename = "./11-monocle/13-monocle_cell_trajectory_for_COL11A1.pdf", width = 4.8, height = 3.5)


plot_cell_trajectory(my_cds_subset, color_by = "seurat_clusters")

BEAM_branch1 <- BEAM(cds = cds, branch_point = 1, cores = 4)



#### 以下代码由于版本问题，无法运行----
my_cds_subset = cds
if(file.exists('BEAM_res.Rdata')){
  
  load(file = 'BEAM_res.Rdata')
}else{
  
  # 这个步骤超级耗费时间
  # 不知道为什么在Mac电脑无法开启并行计算了 
  # 不过我测试了在Windows 电脑设置cores = 4是可以的
  # 如果你是Mac电脑，自己修改 cores = 1 即可 
  BEAM_branch1 <- BEAM(my_cds_subset, branch_point = 1, cores = 4)
  BEAM_branch1 <- BEAM_branch1[order(BEAM_branch1$qval),]
  BEAM_branch1 <- BEAM_branch1[,c("gene_short_name", "pval", "qval")]
  head(BEAM_branch1) 
  
  BEAM_branch2 <- BEAM(my_cds_subset, branch_point = 2, cores = 4)
  BEAM_branch2 <- BEAM_branch2[order(BEAM_branch2$qval),]
  BEAM_branch2 <- BEAM_branch2[,c("gene_short_name", "pval", "qval")]
  head(BEAM_branch2)
  
  save(BEAM_branch1,BEAM_branch2,file = 'BEAM_res.Rdata')
  
}


# 使用全部的基因进行绘图 
BEAM_res = BEAM_branch1
my_branched_heatmap <- plot_genes_branched_heatmap(
  my_cds_subset[row.names(subset(BEAM_res, qval < 1e-4)),],
  branch_point = 1,
  num_clusters = 4, 
  use_gene_short_name = TRUE,
  show_rownames = F,
  return_heatmap = TRUE)

pdf('monocle_BEAM_branch1_heatmap.pdf')
print(my_branched_heatmap$ph)
dev.off()

BEAM_res = BEAM_branch2
my_branched_heatmap <- plot_genes_branched_heatmap(
  my_cds_subset[row.names(subset(BEAM_res, qval < 1e-4)),],
  branch_point = 1,
  num_clusters = 4, 
  use_gene_short_name = TRUE,
  show_rownames = F,
  return_heatmap = TRUE)

pdf('monocle_BEAM_branch2_heatmap.pdf')
print(my_branched_heatmap$ph)
dev.off()



head(my_branched_heatmap$annotation_row)
table(my_branched_heatmap$annotation_row$Cluster) 
my_row <- my_branched_heatmap$annotation_row
my_row <- data.frame(cluster = my_row$Cluster,
                     gene = row.names(my_row),
                     stringsAsFactors = FALSE)

head(my_row[my_row$cluster == 3,'gene']) 

my_gene <- row.names(subset(fData(my_cds_subset),
                            gene_short_name %in% head(my_row[my_row$cluster == 1,'gene'])))
my_gene
# plot genes that are expressed in a branch dependent manner
plot_genes_branched_pseudotime(my_cds_subset[my_gene,],
                               branch_point = 1,
                               ncol = 1)

plot_genes_branched_pseudotime(my_cds_subset[my_gene,],
                               branch_point = 2,
                               ncol = 1)
# 后面的批量绘图，意义不大 
names(pData(my_cds_subset))
head(pData(my_cds_subset))

plot_genes_jitter(my_cds_subset[gene_to_cluster,],
                  grouping = "Cluster",
                  color_by = "Cluster",
                  nrow=  10,
                  ncol = NULL )

ggsave('monocle_top50_subCluster.pdf',height = 42)
plot_genes_in_pseudotime(my_cds_subset[head(gene_to_cluster,25),])
ggsave('monocle_top50_pseudotime.pdf',height = 49)


write.csv(my_pseudotime_de,file = 'my_pseudotime_de.csv')















































































































## 各维度做拟时序分析
p1=plot_cell_trajectory(cds, color_by = "Cluster")  + scale_color_nejm() 
p1
ggsave(p1, filename = "./6-monocle/trajectory_by_cluster.pdf", width = 4.5, height = 3.5)

p2=plot_cell_trajectory(cds, color_by = "Pseudotime")  
p2
ggsave(p2, filename = "./6-monocle/trajectory_by_Pseudotime.pdf", width = 4.5, height = 3.5)

p3=plot_cell_trajectory(cds, color_by = "orig.ident")  + scale_color_nejm()
p3
ggsave(p3, filename = "./6-monocle/trajectory_by_State.pdf", width = 4.5, height = 3.5)

p4=plot_cell_trajectory(cds, color_by = "Cluster")  + facet_wrap("~orig.ident", nrow = 1) + scale_color_nejm()
p4
ggsave(p4, filename = "./6-monocle/trajectory_by_cluster_State.pdf", width = 14, height = 3)

## IGF2BP2做拟时序分析
p5=plot_cell_trajectory(cds, color_by = "Cluster", markers = "MMP11")  + facet_wrap("~orig.ident", nrow = 1) + scale_color_nejm()
p5
ggsave(p5, filename = "./6-monocle/trajectory_by_Igf2bp2_cluster_State.pdf", width = 14, height = 3)

p6=plot_cell_trajectory(cds, color_by = "celltype", markers = "Usp7")  + facet_wrap("~orig.ident", nrow = 1) + scale_color_nejm()
p6
ggsave(p6, filename = "./6-monocle/trajectory_by_Usp7_cluster_State.pdf", width = 14, height = 3)

## IGF2BP2做拟时序分析
p7=plot_cell_trajectory(cds, markers = "Igf2bp2", 
                        use_color_gradient = T) + facet_wrap("~orig.ident", nrow = 1)
p7
ggsave(p7, filename = "./6-monocle/trajectory_by_Igf2bp2_TIME.pdf", width = 14, height = 3)

p8=plot_cell_trajectory(cds, markers = "Usp7", 
                        use_color_gradient = T) + facet_wrap("~orig.ident", nrow = 1)
p8
ggsave(p8, filename = "./6-monocle/trajectory_by_Usp7_TIME.pdf", width = 14, height = 3)





















