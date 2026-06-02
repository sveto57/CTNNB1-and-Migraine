rm(list = ls())
gc()
options(bitmapType='cairo')
setwd("/s1/wuqing/PBS_IS/outputs/IS_naive.6.24")
.libPaths("/home/wuqing/anaconda3/envs/R_4.4/lib/R/library")
library(Seurat)
library(dplyr)
library(tidyverse)
library(Matrix) 
library(data.table) 
library(patchwork)
library(cowplot)
library(ggplot2)
seurat <- readRDS("/s1/wuqing/PBS_IS/outputs/IS_naive.6.24/5.seurat_annotation1.rds")
DimPlot(seurat,reduction = "umap",group.by = "annotation1")
table(seurat$annotation1)
markers <- c("Mgp",#fibroblast_c1
             "Dcn",#fibroblast_c2
             #"",#fibroblast_c3
             "Rbfox3",#neuron
             "Mog",#oligodendrocyte
             "Fabp7","Apoe",#Satglia_C1
             "Ntsr2",#Satglia_C2
             "Mpz",#Schwann 
             "Igfbp7"#vascular
             )
                       
DotPlot(seurat_obj,features =markers ,group.by = "annotation1")
FeaturePlot(seurat_obj,c("Rbfox3","Apoe","Ntsr2","Mog"))
