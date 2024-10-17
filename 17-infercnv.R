
# 设置环境
rm(list = ls())
options(stringsAsFactors = F)
library(Seurat)
library(tidyverse)
library(ggplot2)
library(infercnv)
library(ComplexHeatmap)
library(ggpubr)
library(AnnoProbe)
library(future)
library(ggsci)

# 制作参考CNV集
load(file = "./4-celltype/sce.har.db.Rdata")
Idents(sce.har.db) = sce.har.db@meta.data$celltype

## T细胞参考
Tref = subset(x = sce.har.db, idents = "T cell")
Tref = subset(x = Tref, cells = sample(colnames(Tref), 800))
TrefMat = as.data.frame(GetAssayData(Tref, slot = "counts", assay = "RNA"))

## Endocyte参考
Endoref = subset(x = sce.har.db, idents = "Endothelial")
Endoref = subset(x = Endoref, cells = sample(colnames(Endoref), 800))
EndorefMat = as.data.frame(GetAssayData(Endoref, slot = "counts", assay = "RNA"))

## 保存参考集
save(TrefMat, EndorefMat, file = "./17-infercnv/reference_mat.Rdata") 


# 制作表达矩阵
rm(list = ls())
options(stringsAsFactors = F)
load(file = "./4-celltype/sce.har.db.Rdata")
load(file = "./17-infercnv/reference_mat.Rdata")
Idents(object = sce.har.db) = sce.har.db@meta.data$celltype
sce.har.db = subset(x = sce.har.db, downsample = 1000)
sce.har.db = subset(sce.har.db, idents = c("Ductal", "Acinar", "Endocrine", "Fibroblast", "Stellate", "Macrophage", "B cell", "Plasma cell"))

epiMat = as.data.frame(GetAssayData(sce.har.db, slot = 'counts', assay = 'RNA'))

sce.har.db$celltype = as.character(sce.har.db$celltype)

dat = cbind(epiMat, TrefMat, EndorefMat)
groupinfo = data.frame(v1 = colnames(dat),
                       v2 = c(sce.har.db$celltype, rep("Tref", 800), rep("Endoref", 800)))

# 制作参考基因
geneInfor = annoGene(IDs = rownames(dat), ID_type = "SYMBOL", species = "human")
colnames(geneInfor)
geneInfor = geneInfor[with(geneInfor, order(chr, start)), c(1, 4:6)]
geneInfor = geneInfor[!duplicated(geneInfor[,1]),]
length(unique(geneInfor[,1]))
head(geneInfor)

## 改变染色体排列
geneInfor$chr = factor(geneInfor$chr, levels = c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY"))
geneInfor = geneInfor[order(geneInfor$chr),]
table(geneInfor$chr)
dat = dat[rownames(dat) %in% geneInfor[,1],]
dat = dat[match(geneInfor[,1], rownames(dat)),]
dim(dat)

## 确认顺序
identical(rownames(dat), geneInfor$SYMBOL)
identical(colnames(dat), groupinfo$v1)

## 保存数据
write.table(dat, file = "./17-infercnv/expFile.txt", sep = '\t', quote = F)
write.table(geneInfor, file = "./17-infercnv/geneFile.txt", sep = '\t', quote = F, col.names = F, row.names = F)
write.table(groupinfo, file = "./17-infercnv/groupFiles.txt", sep = '\t', quote = F, col.names = F, row.names = F)

# 运行infercnv
rm(list = ls())
expFile = "./17-infercnv/expFile.txt"
groupFiles = "./17-infercnv/groupFiles.txt"
geneFile = "./17-infercnv/geneFile.txt"

infercnv_obj = CreateInfercnvObject(raw_counts_matrix = expFile,
                                    annotations_file = groupFiles,
                                    delim = "\t",
                                    gene_order_file = geneFile,
                                    ref_group_names = c("Tref", "Endoref"))  ## 这个取决于自己的分组信息里面的

plan("multicore", workers = 4)

infercnv_obj2 = infercnv::run(infercnv_obj,
                              cutoff = 0.1, # cutoff=1 works well for Smart-seq2, and cutoff=0.1 works well for 10x Genomics
                              out_dir = "./17-infercnv/out",  # dir is auto-created for storing outputs
                              cluster_by_groups = T ,   # cluster
                              hclust_method = "ward.D2", 
                              plot_steps = F,
                              denoise = T,
                              HMM = F)

save(infercnv_obj2, file = "./17-infercnv/infercnv_obj2.Rdata")

infercnv::plot_cnv(infercnv_obj2, #上两步得到的infercnv对象
                   output_filename = "./17-infercnv/better_plot",
                   cluster_by_groups = T,
                   output_format = "pdf",
                   custom_color_pal =  color.palette(c("blue","white","red"), c(2, 2)))


# 自定义可视化

## 第一步：提取inferCNV的结果   
obs = infercnv_obj2@expr.data

## 个obs是每一个基因在每一个细胞的拷贝信息,相当于该基因的拷贝量
## 可以通过定义obs的元素数值来让差异变大，在后面画图就能够更大的差异，也可以不运行，更真实体现拷贝数（但是可能没有啥差异）,根据最大最小值来定义
max(obs)
min(obs)
obs[obs < 0.85] <- 1   #把0.6-0.7定义为数值2  后面依此类推
obs[obs >= 0.85 & obs <= 1.05] <- 0
obs[obs > 1.05] <- 1

## 把obs的每一个基因拷贝量加起来，就是这个细胞的总拷贝数obs
score = as.data.frame(colSums(obs))

## 提取meta信息
meta = read.table(file = "./17-infercnv/groupFiles.txt", sep = '\t')  #提取该细胞的其他的meta信息

## 将meta信息添加给score
score = rownames_to_column(score)
score = merge(score, meta, by.x = "rowname", by.y = "V1")    #这里会可能损失一些细胞
score$V2 = factor(score$V2, levels = c("Ductal", "Acinar", "Endocrine", "Fibroblast", "Stellate", "Macrophage", "B cell", "Plasma cell", "Tref", "Endoref"))

## 可视化1-小提琴图
ggplot(data = score, mapping = aes(x = V2, y = colSums(obs), fill = V2))+
  geom_violin()+
  geom_boxplot()+
  theme_bw()
ggsave(filename = "./17-infercnv/1-cnvscore.pdf", width = 5.2, height = 2.7)

## 可视化2-热图
#第一步
#我们先解决细胞注释的表格
#由于我们要定义CNV的程度，所以首先要将数值变成注释
#加change列,标记inferCNV的变化
max(score$`colSums(obs)`)  # 根据最大最小值来定义
min(score$`colSums(obs)`)  # 一般变化倍数大于2倍的可以定义为拷贝数变异
k1 = score$`colSums(obs)` < 1851*2; table(k1)      #这里的250是我min(score$`colSums(obs)`)的结果 
k2 = score$`colSums(obs)` >= 1851*2.5; table(k2)
score <- mutate(score,change = ifelse(k1,"normal",ifelse(k2,"High","Low")))  #要看懂这段代码  ifelse的条件句  假如符合K1标准为normal，不符合就按照ifelse(k2,"High","Low")这个条件再分，假如符合K2标准就High，不符合就Low
table(score$change)   #这样就做好了一个anno

rownames(score) <- score$rowname
score <- score[order(score$`colSums(obs)`),]   #对数据框进行升序排列  如果加decreasing = T就实现降序排列,这样就会对后面热图的细胞顺序进行排序，为什么想对细胞排序呢，因为我想让CNV属于normal，High，Low的细胞排一起，
score <- score[,-c(1,2)]  #删除多余的

#第二步，解决细胞每一个基因的表达量
#由于前面score跟meta merge了，所以这里要对dat进行匹配
dat = infercnv_obj2@count.data
dat = dat[,match(rownames(score), colnames(dat))]
identical(rownames(score), colnames(dat))#这样每一个细胞的anno和它的基因表达量就对应了，记住一定要对应，要想明白哈，有些东西虽然不报错，但是你对应错了，结果也是错的。
dat <- log2(dat+1)#这一步的目的是给表达矩阵取log2，增加基因之间的差异

#第三步 画图啦
#想自定义annotation的颜色       
ann_colors = list(
  change = c(normal="black", Low="white",High="red"))   #这里也可以改变多个anno的颜色  

#画热图
library(pheatmap)
pheatmap(dat,
         show_colnames =F,
         show_rownames = F,
         scale = "row",
         cluster_cols = F,                                        #这些参数建议去Google
         cluster_rows = F, 
         annotation_col=score,
         annotation_colors = ann_colors,
         colorRampPalette(c("#00FF00", "white", "#DC143C"))(75),   #改热图的颜色
         breaks = seq(-3,3,length.out = 100))
































