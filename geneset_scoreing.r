#Drd基因集打分
#Drd基因集打分
library(Seurat)
library(AUCell)
library(msigdbr)
library(dplyr)
library(ggplot2)
library(clusterProfiler)
library(org.Mm.eg.db)
library(GSEABase)
library(GSVA)
rm(list = ls())
gc()
seurat <- readRDS("/s1/wuqing/PBS_IS/outputs/IS_naive.6.24/5.seurat_annotation1.rds")
table(seurat$annotation1,seurat$time)
DimPlot(seurat,group.by = "annotation1")

seurat_neuronal <- subset(seurat, subset = seurat_clusters %in% c("1","2","3","5","6","10","12","13","14"))
seurat_neuronal <- RunPCA(seurat_neuronal, features = VariableFeatures(object = seurat_neuronal) )
ElbowPlot(seurat_neuronal,ndims = 50)
seurat_neuronal <- RunUMAP(seurat_neuronal, reduction = "pca", dims = 1:20)
seurat_neuronal <- FindNeighbors(seurat_neuronal,reduction = "pca", dims = 1:20)
seurat_neuronal <- FindClusters(seurat_neuronal, resolution = 0.5)#12簇
DimPlot(seurat_neuronal,label = T)
Drd <- list(c("Drd1","Drd2","Drd3","Drd4","Drd5",
              "Cacng4","Cdk5r2","Efnb3","Kif1a","Map6",
              "Phf21b","Riiad1","Tagln3","Vangl2"))
migraine <- list(c("Ace","Atp1a2","Cacna1a","Calca","Cfap58",
                   "Mthfr","Notch3","Prrt3","Scn1a","Slc6a4"))
Score_seurat_neuronal <- AddModuleScore(seurat,
                                        features = c(Drd,migraine),
                                        ctrl = 100,      
                                        name = c("Drd","migraine"))   
colnames(Score_seurat_neuronal@meta.data)
VlnPlot(Score_seurat_neuronal,features = c("Drd1","migraine2"),pt.size = 0,adjust = 2,
        split.by = "orig.ident")
DotPlot(Score_seurat_neuronal, features = c("Drd1","migraine2"),split.by = "orig.ident")
#umap图可视化评分结果
library(ggplot2)
mydata <- FetchData(Score_seurat_neuronal,vars = c("umap_1","umap_2","migraine2","annotation1"))
a <- ggplot(mydata,aes(x=umap_1,y=umap_2,colour=migraine2))+
  geom_point(size=1)+scale_color_gradientn(values = seq(0,1,0.2),
                                           colours = c("blue","grey","white","#CC3333"))
a+theme_bw()+theme(panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank(),axis.line = element_line(colour = "black"),
                   panel.border = element_rect(fill = NA,colour = "black",size = 1,linetype = "solid"))

