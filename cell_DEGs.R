#单细胞分组间差异分析Group-wise Differential Expression Analysis
#比较同一细胞类型内部比较不同条件（疾病健康、实验对照）的细胞所表现的差异基因表达模式，基于seurat5
#大类细胞####
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
colors <- c("#BFDFD2","#51999F","#4198AC","#7BC0CD","#DBCB92","#ECB66C","#EA9E58","#ED8D5A","#F1837A")
seurat <- readRDS("/s1/wuqing/PBS_IS/outputs/IS_naive.6.24/5.seurat_annotation1.rds")
DimPlot(seurat,reduction = "umap",group.by = "annotation1",label = T,label.box = T,
        cols = colors,label.size = 2)
table(seurat$orig.ident)
table(seurat$time)
DimPlot(seurat,reduction = "umap",group.by = "annotation1",label = T,label.box = T,
        split.by = "orig.ident",cols = colors,label.size = 2)
seurat$group <- seurat$orig.ident
table(seurat$group,seurat$annotation1)#各组细胞不能为0
Idents(seurat) <- "annotation1"
type <- unique(seurat$annotation1)



#循环：在每一类细胞中的分组间进行差异分析
rm(list = ls())
gc()
setwd("/s1/wuqing/PBS_IS/outputs/IS_naive.6.24")
rownames(r.deg) <- r.deg$X
dir.create("01_组间差异分析_IS", showWarnings = FALSE, recursive = TRUE)
r.deg <- data.frame()
for (i in 1:length(type)) {
  deg <- FindMarkers(seurat,
                     ident.1 = "IS",#处理组"IS_24"
                     ident.2 = "Naive",#对照组"naive"
                     group.by = "group",
                     subset.ident = type[i],
                     min.pct=0.25)
  write.csv(deg,file = paste0("01_组间差异分析_IS/",type[i],"deg.csv"))
  deg$celltype = type[i]
  deg$unm = i-1
  r.deg = rbind(deg,r.deg)}
table(r.deg$celltype)
colnames(r.deg)[which(colnames(r.deg) == "celltype")] <- "cluster"
write.csv(r.deg,file = "/s1/wuqing/PBS_IS/outputs/IS_naive.6.24/01_组间差异分析_IS/celltype_DEGs.csv")
bulk.DEG <- read.csv("/s1/wuqing/Bulk.RNA/allDEGs_p0.05.txt",sep = "")
sc.deg <-subset(r.deg,r.deg$p_val_adj<0.05) 
a <- sc.deg[intersect(rownames(sc.deg),rownames(bulk.DEG)),]
b <- bulk.DEG[intersect(rownames(sc.deg),rownames(bulk.DEG)),]
write.csv(a,file = "/s1/wuqing/PBS_IS/outputs/IS_naive.6.24/01_组间差异分析_IS/bulk.sc.intersection_DEGs.sc.csv")
write.csv(b,file = "/s1/wuqing/PBS_IS/outputs/IS_naive.6.24/01_组间差异分析_IS/bulk.sc.intersection_DEGs.bulk.csv")




#绘图前准备
rm(list = ls())
gc()
setwd("/s1/wuqing/PBS_IS/outputs/IS_naive.6.24")
r.deg <- read.table("./01_组间差异分析_IS/celltype_DEGs.csv",
                  row.names = 1,header=T,sep = ",")
rownames(r.deg) <- r.deg$X
r.deg <- r.deg[,-1]
s.deg <- subset(r.deg,p_val_adj < 0.05&abs(avg_log2FC)>0.5)
table(s.deg$cluster)
s.deg$threshold <- as.factor(ifelse(s.deg$avg_log2FC>1,"up","down"))
table(s.deg$threshold)
dim(s.deg)
s.deg$adj_p_signi <- as.factor(ifelse(s.deg$p_val_adj<0.01,"Highly","Lowly"))
table(s.deg$adj_p_signi)
s.deg$thr_signi <- paste0(s.deg$threshold,"_",s.deg$adj_p_signi)
write.csv(s.deg,"01_组间差异分析_IS/deg_all_IS.Naive.csv")


#火山图可视化
top_up_label <- s.deg%>%
  subset(.,threshold%in%"up")%>%
  group_by(cluster)%>%
  top_n(n=5,wt=avg_log2FC)%>%
  as.data.frame()
top_down_label <- s.deg%>%
  subset(.,threshold%in%"down")%>%
  group_by(cluster)%>%
  top_n(n=-5,wt=avg_log2FC)%>%
  as.data.frame()
top_label <- rbind(top_up_label,top_down_label)
library(scRNAtoolVis)
library(ggrepel)
library(dplyr)
colnames(r.deg)
head(r.deg)
r.deg$gene <- rownames(r.deg)
jjVolcano(diffData=r.deg,
          tile.col=colors[1:9],#细胞类型数量
          pSize=0.4,
          legend.position=c(0.1,0.9),
          celltypeSize=2,
          topGeneN=5)+
  labs(title = "IS.Naive")

p <- scRNAtoolVis::markerVocalno(markers = r.deg,
                                 topn = 5,
                                 labelCol = colors)
setwd("/s1/wuqing/PBS_IS/outputs/IS_naive.6.24/01_组间差异分析_IS")
ggsave("markerVolcano.pdf", plot = p, width = 10, height = 8)














#neuron_sub细胞####
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
colors <- c("#BFDFD2","#51999F","#4198AC","#7BC0CD","#DBCB92","#ECB66C","#EA9E58","#ED8D5A","#F1837A")
seurat <- readRDS("/s1/wuqing/PBS_IS/outputs/IS_naive.6.24/5.neuron_subtype.rds")
DimPlot(seurat,reduction = "umap",group.by = "subtype",label = T,label.box = T,
        cols = colors,label.size = 2)
table(seurat$orig.ident)
table(seurat$time)
DimPlot(seurat,reduction = "umap",group.by = "subtype",label = T,label.box = T,
        split.by = "orig.ident",cols = colors,label.size = 2)
seurat$group <- seurat$orig.ident
table(seurat$group,seurat$subtype)#各组细胞不能为0
Idents(seurat) <- "subtype"
type <- unique(seurat$subtype)



#循环：在每一类细胞中的分组间进行差异分析
rm(list = ls())
gc()
setwd("/s1/wuqing/PBS_IS/outputs/IS_naive.6.24/01_组间差异分析_IS/neuron_subtype_DEGs")
dir.create("neuron_subtype_DEGs", showWarnings = FALSE, recursive = TRUE)
r.deg <- data.frame()
for (i in 1:length(type)) {
  deg <- FindMarkers(seurat,
                     ident.1 = "IS",#处理组"IS_24"
                     ident.2 = "Naive",#对照组"naive"
                     group.by = "group",
                     subset.ident = type[i],
                     min.pct=0.25)
  write.csv(deg,file = paste0("neuron_subtype_DEGs",type[i],"deg.csv"))
  deg$celltype = type[i]
  deg$unm = i-1
  r.deg = rbind(deg,r.deg)}
table(r.deg$celltype)
r.deg$X <- rownames(r.deg) 
colnames(r.deg)[which(colnames(r.deg) == "celltype")] <- "cluster"
setwd("/s1/wuqing/PBS_IS/outputs/IS_naive.6.24/01_组间差异分析_IS/neuron_subtype_DEGs")
write.csv(r.deg,file = "all_subtype_DEGs.csv")
sc.deg <-subset(r.deg,r.deg$p_val_adj<0.05) 
write.csv(sc.deg,file = "all_subtype_DEGs_0.05p.csv")


#绘图前准备
rm(list = ls())
gc()
setwd("/s1/wuqing/PBS_IS/outputs/IS_naive.6.24/01_组间差异分析_IS/neuron_subtype_DEGs")
r.deg <- read.table("all_subtype_DEGs.csv",
                    row.names = 1,header=T,sep = ",")
rownames(r.deg) <- r.deg$X
r.deg$X <- NULL
p.deg <- subset(r.deg,p_val_adj < 0.05)
s.deg <- subset(r.deg,p_val_adj < 0.05&abs(avg_log2FC)>0.5)
table(p.deg$cluster)
table(s.deg$cluster)
s.deg$threshold <- as.factor(ifelse(s.deg$avg_log2FC>0,"up","down"))
table(s.deg$threshold)
dim(s.deg)
s.deg$adj_p_signi <- as.factor(ifelse(s.deg$p_val_adj<0.01,"Highly","Lowly"))
table(s.deg$adj_p_signi)
s.deg$thr_signi <- paste0(s.deg$threshold,"_",s.deg$adj_p_signi)
intersection <- c("Pcmtd2","R3hdm1","Bag3","Cacna2d1","Nat8l","Ttc3","Sfswap","Arfgef3","Htt","Pcdh7","Gria2","Hdac11","Tbc1d16","Gfra1","Ndufa6","Mei4","Srsf11","Rbm6","Ebf3" ,"Lgals1","Tmod3")
intersection_gene <- r.deg[intersection,]
intersection_gene <- intersection_gene[intersection_gene$p_val_adj<0.05,]
#火山图可视化
top_up_label <- s.deg%>%
  subset(.,threshold%in%"up")%>%
  group_by(cluster)%>%
  top_n(n=3,wt=avg_log2FC)%>%
  as.data.frame()
top_down_label <- s.deg%>%
  subset(.,threshold%in%"down")%>%
  group_by(cluster)%>%
  top_n(n=-3,wt=avg_log2FC)%>%
  as.data.frame()
top_label <- rbind(top_up_label,top_down_label)
library(scRNAtoolVis)
library(ggrepel)
library(dplyr)
colnames(r.deg)
head(r.deg)
r.deg$gene <- rownames(r.deg)
colors <- c("#BFDFD2","#51999F","#4198AC","#7BC0CD","#DBCB92","#ECB66C","#EA9E58","#ED8D5A","#F1837A")
p.deg <- r.deg[r.deg$p_val_adj<0.05,]
jjVolcano(diffData=p.deg,
          tile.col=colors[1:7],#细胞类型数量
          pSize=0.4,
          legend.position=c(0.1,0.9),
          celltypeSize=1)+
  labs(title = "neuron_subtype")
p <- scRNAtoolVis::markerVocalno(markers = r.deg,
                                 topn = 3,
                                 labelCol = colors)
p
setwd("/s1/wuqing/PBS_IS/outputs/IS_naive.6.24/01_组间差异分析_IS/neuron_subtype_DEGs")
ggsave("markerVolcano.pdf", plot = p, width = 10, height = 8)


#富集####
deg_NP.Mrgprd <- r.deg[r.deg$cluster=="NP_Mrgprd+",]
write.csv(deg_NP.Mrgprd,"neuron_subtype_DEGsNP_Mrgprd+deg.csv")
deg_NP.Mrgprd <- deg_NP.Mrgprd[deg_NP.Mrgprd$p_val_adj< 0.05,]
up.deg_NP.Mrgprd <- deg_NP.Mrgprd[deg_NP.Mrgprd$avg_log2FC> 0.5,]
down.deg_NP.Mrgprd <- deg_NP.Mrgprd[deg_NP.Mrgprd$avg_log2FC< -0.5,]
library(tidyverse)
library(clusterProfiler)
library(enrichplot)
library(org.Mm.eg.db)
library(dplyr)
library(msigdbr)
up.gene_ids <- bitr(up.deg_NP.Mrgprd$X, 
                 fromType = "SYMBOL", 
                 toType = "ENTREZID", 
                 OrgDb = org.Mm.eg.db)
down.gene_ids <- bitr(down.deg_NP.Mrgprd$X, 
                    fromType = "SYMBOL", 
                    toType = "ENTREZID", 
                    OrgDb = org.Mm.eg.db)
#GO
go_result.up <- enrichGO(gene = up.gene_ids$ENTREZID,
                      OrgDb = org.Mm.eg.db,
                      ont = "ALL",
                      pAdjustMethod = "BH",
                      pvalueCutoff = 0.05,
                      qvalueCutoff = 0.2,
                      readable = TRUE )
setwd("/s1/wuqing/PBS_IS/outputs/IS_naive.6.24/01_组间差异分析_IS/neuron_subtype_DEGs")
saveRDS(go_result.up,"go_result_NP.Mrgprd.up.rds")
Go_result.up <- as.data.frame(go_result.up)
write.csv(Go_result.up,"Go_result_NP.Mrgprd.up.csv")
go_result.down <- enrichGO(gene = down.gene_ids$ENTREZID,
                         OrgDb = org.Mm.eg.db,
                         ont = "ALL",
                         pAdjustMethod = "BH",
                         pvalueCutoff = 0.05,
                         qvalueCutoff = 0.2,
                         readable = TRUE )
saveRDS(go_result.down,"go_result_NP.Mrgprd.down.rds")
Go_result.down <- as.data.frame(go_result.down)
write.csv(Go_result.down,"Go_result_NP.Mrgprd.down.csv")

dotplot(go_result, showCategory = 10, split = "ONTOLOGY") +
  facet_grid(ONTOLOGY ~ ., scales = "free")
barplot(go_result, showCategory = 10)
#KEGG
kegg_result_up <- enrichKEGG(gene = up.gene_ids$ENTREZID,organism = 'mmu',pvalueCutoff = 0.05)
kegg_result_down <- enrichKEGG(gene = down.gene_ids$ENTREZID,organism = 'mmu',pvalueCutoff = 0.05)
KEGG_result_up <- as.data.frame(kegg_result_up)
KEGG_result_down <- as.data.frame(kegg_result_down)
write.csv(KEGG_result_up, "KEGG_result_NP.Mrgprd.up.csv")
write.csv(KEGG_result_down, "KEGG_result_NP.Mrgprd.down.csv")
saveRDS(kegg_result_up,"kegg_result_NP.Mrgprd.up.rds")
saveRDS(kegg_result_down,"kegg_result_NP.Mrgprd.down.rds")

#结果可视化
Go_result.down
Go_result.up
KEGG_result_down
KEGG_result_up
library(clusterProfiler)
dotplot(Go_result.down, showCategory = 15)
clusterProfiler::dotplot(KEGG_result_up, showCategory = 15)
clusterProfiler::barplot() 

#指定通路高亮，横坐标为-log10P
library(ggplot2)
library(dplyr)
setwd("/s1/wuqing/PBS_IS/outputs/IS_naive.6.24/01_组间差异分析_IS/neuron_subtype_DEGs")
#KEGG_result_up####
pathways_keep <- c("Dopaminergic synapse - Mus musculus (house mouse)", 
                   "Glutamatergic synapse - Mus musculus (house mouse)",
                   "Citrate cycle (TCA cycle) - Mus musculus (house mouse)",
                   "Long-term depression - Mus musculus (house mouse)",
                   "HIF-1 signaling pathway - Mus musculus (house mouse)",
                   "Estrogen signaling pathway - Mus musculus (house mouse)")
plot_data_filtered <-  %>%
  filter(Description %in% pathways_keep) %>%
  mutate(log10_padj = -log10(p.adjust))
ggplot(plot_data_filtered, aes(x = log10_padj, y = reorder(Description, log10_padj))) +
  geom_bar(aes(fill = p.adjust), stat = "identity", width = 0.7) +
  scale_fill_gradient(low = "red", high = "blue", name = "P.adjust") +
  labs(x = "-log10(P.adj)",y=NULL,title = "KEGG Pathways of NP_Mrgprd+ upDEGs") +
  theme_bw() +
  theme(axis.text.y = element_text(size = 11, color = "black"),
        axis.text.x = element_text(size = 10, color = "black"),
        axis.title = element_text(size = 12, face = "bold"),
        plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
        legend.position = "right")+
  scale_x_continuous(limits = c(0, 1.4),     # 设定范围
                     breaks = seq(0,1.4,0.1),  # 刻度位置(每隔1)
                     expand = c(0, 0)   )      # 去除边距)
    
#KEGG_result_down####
pathways_keep <- c("Oxidative phosphorylation - Mus musculus (house mouse)", 
                   "IgSF CAM signaling - Mus musculus (house mouse)")
plot_data_filtered <- KEGG_result_down %>%
  filter(Description %in% pathways_keep) %>%
  mutate(log10_padj = -log10(p.adjust))
range(plot_data_filtered$log10_padj)
ggplot(plot_data_filtered, aes(x = log10_padj, y = reorder(Description, log10_padj))) +
  geom_bar(aes(fill = p.adjust), stat = "identity", width = 0.7) +
  scale_fill_gradient(low = "red", high = "blue", name = "P.adjust") +
  labs(x = "-log10(P.adj)",y=NULL,title = "KEGG Pathways of NP_Mrgprd+ downDEGs") +
  theme_bw() +
  theme(axis.text.y = element_text(size = 11, color = "black"),
        axis.text.x = element_text(size = 10, color = "black"),
        axis.title = element_text(size = 12, face = "bold"),
        plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
        legend.position = "right")+
  scale_x_continuous(limits = c(0, 5.7),     
                     breaks = seq(0,5.7,0.5),  
                     expand = c(0, 0)   )      




#GO_result_up####
a <- Go_result.up[Go_result.up$ONTOLOGY %in%c("CC","BP"),] %>%
  filter(grepl("mitochond|mitochondrial|lipid|phosphorylation|complex|fatty acid|metabolic", 
               Description, TRUE))
pathways_keep <- c("negative regulation of Ras protein signal transduction", 
                   "behavioral response to pain",
                  "regulation of calcium ion-dependent exocytosis",
                  "synaptic vesicle endocytosis",
                  "mitochondrion transport along microtubule",
                  "acetyl-CoA biosynthetic process")
plot_data_filtered <- Go_result.up %>%
  filter(Description %in% pathways_keep) %>%
  mutate(log10_padj = -log10(p.adjust))
range(plot_data_filtered$log10_padj)
ggplot(plot_data_filtered, aes(x = log10_padj, y = reorder(Description, log10_padj))) +
  geom_bar(aes(fill = p.adjust), stat = "identity", width = 0.7) +
  scale_fill_gradient(low = "red", high = "blue", name = "P.adjust") +
  labs(x = "-log10(P.adj)",y=NULL,title = "GO Pathways of NP_Mrgprd+ upDEGs") +
  theme_bw() +
  theme(axis.text.y = element_text(size = 11, color = "black"),
        axis.text.x = element_text(size = 10, color = "black"),
        axis.title = element_text(size = 12, face = "bold"),
        plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
        legend.position = "right")+
  scale_x_continuous(limits = c(0, 2.5),     
                     breaks = seq(0,2.5,0.5),  
                     expand = c(0, 0)   )      



#GO_result_down####
a <- Go_result.down%>%filter(
  grepl("mitochond|mitochondrial|lipid|phosphorylation|complex|fatty acid|metabolic", 
  Description, TRUE))
pathways_keep <- c("mitochondrial electron transport, cytochrome c to oxygen", 
                   "mitochondrial ATP synthesis coupled electron transport",
                   "positive regulation of dephosphorylation",
                   "beta-catenin destruction complex",
                   "ATPase regulator activity")
plot_data_filtered <- Go_result.down %>%
  filter(Description %in% pathways_keep) %>%
  mutate(log10_padj = -log10(p.adjust))
range(plot_data_filtered$log10_padj)
ggplot(plot_data_filtered, aes(x = log10_padj, y = reorder(Description, log10_padj))) +
  geom_bar(aes(fill = p.adjust), stat = "identity", width = 0.7) +
  scale_fill_gradient(low = "red", high = "blue", name = "P.adjust") +
  labs(x = "-log10(P.adj)",y=NULL,title = "GO Pathways of NP_Mrgprd+ downDEGs") +
  theme_bw() +
  theme(axis.text.y = element_text(size = 11, color = "black"),
        axis.text.x = element_text(size = 10, color = "black"),
        axis.title = element_text(size = 12, face = "bold"),
        plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
        legend.position = "right")+
  scale_x_continuous(limits = c(0, 3.9),     
                     breaks = seq(0,3.9,0.5),  
                     expand = c(0, 0)   )   






































