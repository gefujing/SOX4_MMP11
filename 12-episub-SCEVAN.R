
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
library(ggrepel)
library(enrichplot)
library(SCEVAN)
load("./6-episub/episub.Rdata")
load(file = "./0-Rawdata/my_color.Rdata")

# 准备输入数据

# episub-cluster0
episub0 = subset(episub, seurat_clusters == 0)
episub0 = subset(episub0, downsample = 3000)
count_mtx = episub0@assays$RNA@counts
results = pipelineCNA(count_mtx, par_cores = 10, SUBCLONES = F)

# episub-cluster1
episub1 = subset(episub, seurat_clusters == 1)
episub1 = subset(episub1, downsample = 3000)
count_mtx = episub1@assays$RNA@counts
results = pipelineCNA(count_mtx, par_cores = 10, SUBCLONES = F)










































