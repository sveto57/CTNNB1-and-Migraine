rm(list = ls())
gc()
options(bitmapType='cairo')
setwd("/s1/wuqing/PBS_IS/outputs")
.libPaths("/home/wuqing/anaconda3/envs/R_4.4/lib/R/library")
library(Seurat)
library(dplyr)
library(tidyverse)
library(Matrix) 
library(data.table) 
seurat<- readRDS("/s1/wuqing/PBS_IS/outputs/IS_chronological/0.IS.Naive_indrops.rds")
set.seed(123)
naive_samples <- c("Naive_male_rep1", "Naive_male_rep2", "Naive_male_rep3", 
                   "Naive_male_rep4", "Naive_male_rep5", "Naive_male_rep6", 
                   "Naive_male_rep7", "Naive_male_rep8", "Naive_male_rep9")
naive_samples <- sample(naive_samples,2,replace = F)
seurat_obj <- subset(seurat, subset = sample %in% c("IS_24h_male_rep1","IS_24h_male_rep2",
                                                        "IS_6h_male_rep1","IS_6h_male_rep2",
                                                        "Naive_male_rep5","Naive_male_rep4"))
seurat_obj$time <- "naive"
seurat_obj$time <-  ifelse(seurat_obj$sample %in% c("IS_24h_male_rep1","IS_24h_male_rep2"),"IS_24",
                           ifelse(seurat_obj$sample %in% c("IS_6h_male_rep1","IS_6h_male_rep2"),"IS_6","naive"))
table(seurat_obj$time)
saveRDS(seurat_obj,file = "/s1/wuqing/PBS_IS/outputs/IS_naive.6.24/0.IS_6.24_Naive.rds")






#QC####
rm(list = ls())
gc()
options(bitmapType='cairo')
setwd("/s1/wuqing/PBS_IS/outputs")
.libPaths("/home/wuqing/anaconda3/envs/R_4.4/lib/R/library")
library(Seurat)
library(dplyr)
library(ggplot2)
library(tidyverse)
library(Matrix) 
library(data.table) 
library(harmony)  
library(DoubletFinder)
seurat_obj <- read_rds("/s1/wuqing/PBS_IS/outputs/IS_naive.6.24/0.IS_6.24_Naive.rds")
seurat_obj
table(seurat_obj$time)



seurat_obj [["percent.mt"]] <- PercentageFeatureSet(seurat_obj , pattern = "^mt-")
seurat_obj [["percent.rbc"]] <- PercentageFeatureSet(seurat_obj , pattern = "^Hb[ab](-|$)")
seurat_obj
VlnPlot(seurat_obj,pt.size=0,ncol = 2,group.by = "time",
        features = c( "nFeature_RNA", "nCount_RNA", "percent.mt","percent.rbc"))
quantile(seurat_obj$percent.rbc,seq(0,1,0.1))
quantile(seurat_obj$percent.mt,seq(0.9,1,0.01))
seurat_obj#7214个细胞
seurat_obj  <- subset(seurat_obj, 
                          subset = nFeature_RNA>400 & nCount_RNA<40000 & percent.mt<5 & percent.rbc < 1)
seurat_obj#7026个细胞
saveRDS(seurat_obj,file = "/s1/wuqing/PBS_IS/outputs/IS_naive.6.24/1.IS_6.24_Naive_QC.rds")



seurat_obj
seurat_obj <- NormalizeData(seurat_obj, normalization.method = "LogNormalize", scale.factor = 10000)
seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000)
seurat_obj <- ScaleData(seurat_obj)
seurat_obj <- RunPCA(seurat_obj, features = VariableFeatures(object = seurat_obj) )
ElbowPlot(seurat_obj,ndims = 50)
seurat_obj <- RunUMAP(seurat_obj, reduction = "pca", dims = 1:20)
seurat_obj <- FindNeighbors(seurat_obj,reduction = "pca", dims = 1:20)
seurat_obj <- FindClusters(seurat_obj)
seurat_obj


#细胞周期评分
capitalize_genes <- function(genes) {
  return(sapply(genes, function(x) {
    paste0(toupper(substring(x, 1, 1)), tolower(substring(x, 2)))
  }))
}
s.genes <- capitalize_genes(cc.genes$s.genes) 
g2m.genes <- capitalize_genes(cc.genes$g2m.genes)
head(seurat_obj$gene)
seurat_obj <- CellCycleScoring(seurat_obj,
                               s.features = s.genes,g2m.features = g2m.genes,
                               set.ident = TRUE)
seurat_obj@meta.data[1:5,]
table(seurat_obj$Phase)
DimPlot(seurat_obj,group.by = "Phase")
table(seurat_obj$Phase,seurat_obj$seurat_clusters)
saveRDS(seurat_obj, file = "/s1/wuqing/PBS_IS/outputs/IS_naive.6.24/2.IS_Naive_QC_Phase.rds")



#annotation1####
FeaturePlot(seurat_obj,features = c("Rbfox3","Sparc"),reduction = "umap",label = T)
table(seurat_obj$seurat_clusters)
seurat_obj$neuron_or_not <- recode(seurat_obj$seurat_clusters,
                            "0"="nonneuron",
                            "1"="nonneuron",
                            "2"="neuron",
                            "3"="neuron",
                            "4"="neuron",
                            "5"="neuron",
                            "6"="neuron",
                            "7"="nonneuron",
                            "8"="nonneuron",
                            "9"="nonneuron",
                            "10"="neuron",
                            "11"="nonneuron?",
                            "12"="neuron",
                            "13"="neuron",
                            "14"="neuron",
                            "15"="nonneuron",
                            "16"="nonneuron",
                            "17"="neuron",
                            "18"="nonneuron")
saveRDS(seurat_obj,file = "/s1/wuqing/PBS_IS/outputs/IS_naive.6.24/3.seurat_annotation1.rds")       





#low_quality_cluster
rm(list = ls())
gc()
options(bitmapType='cairo')
setwd("/s1/wuqing/PBS_IS/outputs")
.libPaths("/home/wuqing/anaconda3/envs/R_4.4/lib/R/library")
library(Seurat)
library(dplyr)
library(ggplot2)
library(tidyverse)
library(Matrix) 
library(data.table) 
library(harmony)  
library(DoubletFinder)
seurat_obj <- read_rds("/s1/wuqing/PBS_IS/outputs/IS_naive.6.24/3.seurat_annotation1.rds")
unique(seurat_obj$neuron_or_not)
seurat_neuron <- subset(seurat_obj, subset = neuron_or_not == "neuron")
seurat_nonneuron <- subset(seurat_obj, subset = neuron_or_not %in% c("nonneuron","nonneuron?"))


seurat_neuron#3624个细胞
seurat_neuron <- RunPCA(seurat_neuron, features = VariableFeatures(object = seurat_neuron) )
seurat_neuron <- RunUMAP(seurat_neuron, dims = 1:20)
seurat_neuron <- FindNeighbors(seurat_neuron,  dims = 1:20)
seurat_neuron <- FindClusters(seurat_neuron)
allMarkers_neuron <- FindAllMarkers(seurat_neuron,
                             only.pos = TRUE,          
                             logfc.threshold = 0,      
                             min.pct = 0.1)
mito_genes <- grep("^mt-", rownames(seurat_neuron), value = TRUE, ignore.case = TRUE)
clusters_mito_high <- allMarkers_neuron %>%
  filter(gene %in% mito_genes,
         avg_log2FC > 1,
         p_val_adj < 0.05) %>%
  group_by(cluster) %>%
  summarise(n_mito = n(), .groups = "drop") %>%
  filter(n_mito >= 2) %>%
  pull(cluster)%>%
  as.character()
clusters_mito_high
a <- allMarkers_neuron[allMarkers_neuron$p_val_adj<0.05,]
table(a$cluster)

seurat_nonneuron #3402个细胞
seurat_nonneuron <- RunPCA(seurat_nonneuron, features = VariableFeatures(object = seurat_neuron) )
seurat_nonneuron <- RunUMAP(seurat_nonneuron, dims = 1:20)
seurat_nonneuron <- FindNeighbors(seurat_nonneuron,  dims = 1:20)
seurat_nonneuron <- FindClusters(seurat_nonneuron)
allMarkers_nonneuron <- FindAllMarkers(seurat_nonneuron,
                                    only.pos = TRUE,          
                                    logfc.threshold = 0,      
                                    min.pct = 0.1)
mito_genes <- grep("^mt-", rownames(seurat_nonneuron), value = TRUE, ignore.case = TRUE)
clusters_mito_high <- allMarkers_nonneuron %>%
  filter(gene %in% mito_genes,
         avg_log2FC > 1,
         p_val_adj < 0.05) %>%
  group_by(cluster) %>%
  summarise(n_mito = n(), .groups = "drop") %>%
  filter(n_mito >= 2) %>%
  pull(cluster)%>%
  as.character()
clusters_mito_high


markers_c11 <- allMarkers_nonneuron[allMarkers_nonneuron$cluster==11,]
markers_c11 <- markers_c11[markers_c11$p_val_adj<0.05,]
seurat_11 <- subset(seurat_nonneuron,subset = seurat_clusters == 11)
table(seurat_11$orig.ident)
table(seurat_nonneuron$seurat_clusters)
rownames(seurat_11@meta.data)
condition <- intersect(rownames(seurat_obj@meta.data),rownames(seurat_11@meta.data))
seurat_obj$low_quality <- ifelse(rownames(seurat_obj@meta.data) %in% condition,
                                 TRUE,FALSE)
table(seurat_obj$low_quality,useNA="ifany")
DimPlot(seurat_obj,reduction = "umap",group.by = "low_quality")
table(seurat_obj$low_quality,seurat_obj$annotation1)
saveRDS(seurat_obj,file = "/s1/wuqing/PBS_IS/outputs/IS_naive.6.24/4.1.seurat_lowquality.rds")
seurat <- subset(seurat_obj,subset = low_quality == FALSE)
saveRDS(seurat,file = "/s1/wuqing/PBS_IS/outputs/IS_naive.6.24/4.2.seurat_lowquality_rm.rds")






rm(list = ls())
gc()
options(bitmapType='cairo')
setwd("/s1/wuqing/PBS_IS/outputs")
.libPaths("/home/wuqing/anaconda3/envs/R_4.4/lib/R/library")
library(Seurat)
library(dplyr)
library(ggplot2)
library(tidyverse)
library(Matrix) 
library(data.table) 
library(harmony)  
seurat_obj <- readRDS("/s1/wuqing/PBS_IS/outputs/IS_naive.6.24/4.2.seurat_lowquality_rm.rds")
seurat_obj#6932个细胞
seurat_obj <- NormalizeData(seurat_obj, normalization.method = "LogNormalize", scale.factor = 10000)
seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000)
seurat_obj <- ScaleData(seurat_obj)
seurat_obj <- RunPCA(seurat_obj, features = VariableFeatures(object = seurat_obj) )
ElbowPlot(seurat_obj,ndims = 50)
seurat_obj <- RunUMAP(seurat_obj, reduction = "pca", dims = 1:20)
seurat_obj <- FindNeighbors(seurat_obj,reduction = "pca", dims = 1:20)
seurat_obj <- FindClusters(seurat_obj, resolution = c(0.5,0.6,0.7,0.8, 0.9, 1.0, 1.1, 1.2))


#dims选择
library(dplyr)
library(magrittr)
library(patchwork)
seurat_obj <- readRDS("/s1/wuqing/PBS_IS/outputs/IS_naive.6.24/5.seurat_annotation1.rds")
for (i in c(15,20,25,30,35,40)) {
  seurat <- FindNeighbors(seurat,reduction = "pca", dims = 1:i) |> FindClusters(resolution = 0.5)
  seurat <- RunUMAP(seurat, reduction = "pca", dims = 1:i)
  plot_i <- print(DimPlot(seurat,reduction = "umap",label = T,pt.size = 1,label.size = 5,repel = T)+labs(title = i))
  plotname <- paste("plot_",i,sep="")
  assign(plotname,plot_i)
  print(plot_i)
}
p <- plot_15+plot_20+plot_25+plot_30+plot_35+plot_40+plot_layout(ncol = 3)
setwd("/s1/wuqing/PBS_IS/outputs/IS_naive.6.24")
ggsave(p,filename ="dims.15_40.pdf",width = 15,height = 8)
#resolution选择
library(clustree)
library(patchwork)
set.seed(123)
seurat_obj <- FindNeighbors(seurat,reduction = "pca", dims = 1:20)
seurat<- FindClusters(seurat, resolution = seq(0.2,1.2,0.1))
a <- clustree(seurat,prefix="RNA_snn_res.")+coord_flip()
b <- a+plot_layout(widths = c(3:1))
ggsave(b,filename="clustree.pdf",width=12,height = 9)
colnames(seurat@meta.data)
table(seurat$RNA_snn_res.0.5)

p1 <- DimPlot(seurat,reduction = "umap",group.by = "RNA_snn_res.0.2",
              label = T,repel = F,shuffle = T)
p2 <- DimPlot(seurat,reduction = "umap",group.by = "RNA_snn_res.0.3",
              label = T,repel = F,shuffle = T)
p3 <- DimPlot(seurat,reduction = "umap",group.by = "RNA_snn_res.0.4",
              label = T,repel = F,shuffle = T)
p4 <- DimPlot(seurat,reduction = "umap",group.by = "RNA_snn_res.0.5",
              label = T,repel = F,shuffle = T)
p5 <- DimPlot(seurat,reduction = "umap",group.by = "RNA_snn_res.0.6",
              label = T,repel = F,shuffle = T)
p6 <- DimPlot(seurat,reduction = "umap",group.by = "RNA_snn_res.0.7",
              label = T,repel = F,shuffle = T)
p7 <- DimPlot(seurat,reduction = "umap",group.by = "RNA_snn_res.0.8",
              label = T,repel = F,shuffle = T)
p8 <- DimPlot(seurat,reduction = "umap",group.by = "RNA_snn_res.0.9",
              label = T,repel = F,shuffle = T)
p9 <- DimPlot(seurat,reduction = "umap",group.by = "RNA_snn_res.0.9",
              label = T,repel = F,shuffle = T)
p10 <- DimPlot(seurat,reduction = "umap",group.by = "RNA_snn_res.1",
              label = T,repel = F,shuffle = T)
p11 <- DimPlot(seurat,reduction = "umap",group.by = "RNA_snn_res.1.1",
              label = T,repel = F,shuffle = T)
p12 <- DimPlot(seurat,reduction = "umap",group.by = "RNA_snn_res.1.2",
              label = T,repel = F,shuffle = T)
p <- p1+p2+p3+p4+p5+p6+p7+p8+p9+p10+p11+p12+plot_layout(ncol = 5)
ggsave(p,filename="RNA_snn_res.0.2_1.2.pdf",width = 25,height = 8)





table(seurat_obj$neuron_or_not,seurat_obj$seurat_clusters)
nonneuron_markers <- c("Cd74",#immune "Ptprc",
                       "Igfbp7",#vascular
                       "Dcn","Mgp",#fibroblast
                       "Apoe",#Satglia "Fabp7",
                       "Mpz",#Schwann 
                       "Mog","Ermn")#oligodendrocyte
c("Mpz","Gldn",#Schwann_M
  "Scn7a")#Schwann_N
FeaturePlot(seurat_obj,features=c("Rbfox3","Sparc"),reduction = "umap",label = T)
FeaturePlot(seurat_obj,features=c("Mpz","Gldn","Scn7a"),reduction = "umap",label = T)
DotPlot(seurat_obj,features = nonneuron_markers)+RotatedAxis()
seurat_obj$annotation1 <- recode(seurat_obj$seurat_clusters,
                                 "0"="schwann",
                                 "1"="neuron",
                                 "2"="satglia",
                                 "3"="neuron",
                                 "4"="neuron",
                                 "5"="neuron",
                                 "6"="neuron",
                                 "7"="satglia?",
                                 "8"="vascular",
                                 "9"="fibroblast_c1",
                                 "10"="fibroblast_c2",
                                 "11"="neuron",
                                 "12"="oligodendrocyte",
                                 "13"="neuron",
                                 "14"="neuron",
                                 "15"="neuron",
                                 "16"="satglia",
                                 "17"="fibroblast_c3",
                                 "18"="neuron",
                                 "19"="satglia")
seurat_obj <- FindSubCluster(seurat_obj, cluster = 7,
                             graph.name = "RNA_snn",
                             subcluster.name = "sub.cluster.7")
DimPlot(seurat_obj,group.by = "sub.cluster.7",label = T)
DotPlot(seurat_obj,features = nonneuron_markers,group.by = "sub.cluster.7")+RotatedAxis()
seurat_obj$annotation1 <- recode(seurat_obj$sub.cluster.7,
                                 "0"="schwann",
                                 "1"="neuron",
                                 "2"="satglia_c1",
                                 "3"="neuron",
                                 "4"="neuron",
                                 "5"="neuron",
                                 "6"="neuron",
                                 "7_0"="satglia_c1",
                                 "7_1"="satglia_c1",
                                 "7_2"="schwann",
                                 "8"="vascular",
                                 "9"="fibroblast_c1",
                                 "10"="fibroblast_c2",
                                 "11"="neuron",
                                 "12"="oligodendrocyte",
                                 "13"="neuron",
                                 "14"="neuron",
                                 "15"="neuron",
                                 "16"="satglia_c2",
                                 "17"="fibroblast_c2",
                                 "18"="neuron",
                                 "19"="satglia_c1")
table(seurat_obj$annotation1)
table(seurat_obj$neuron_or_not,seurat_obj$annotation1)
saveRDS(seurat_obj,file = "/s1/wuqing/PBS_IS/outputs/IS_naive.6.24/5.seurat_annotation1.rds")       
DimPlot(seurat_obj,group.by = "annotation1",split.by = "orig.ident")

neuron_markers <-c("Slc17a6",#谷氨酸能神经元
                   "Slc32a1",#GABA神经元
                   "Th",#DA能神经元
                   "Chat")#胆碱能 
DotPlot(seurat_obj,features = neuron_markers)+RotatedAxis()
DEG <- FindAllMarkers(seurat_obj,
                      min.pct=0.25)


#cell proportion####
rm(list = ls())
gc()
seurat <- readRDS("/s1/wuqing/PBS_IS/outputs/IS_naive.6.24/5.seurat_annotation1.rds")
table(seurat$annotation1)
table(seurat$orig.ident,seurat$annotation1)
table(seurat$time,seurat$annotation1)
meta <- FetchData(seurat,vars = c("orig.ident", "annotation1","time","sample"))
meta$annotation1 <- factor(meta$annotation1,
                           levels = c("fibroblast_c1","fibroblast_c2","fibroblast_c3","vascular",
                                      "neuron","oligodendrocyte", "satglia_c1", "satglia_c2","schwann"))
meta$orig.ident<- factor(meta$orig.ident,levels = c("Naive","IS"))
meta$time<- factor(meta$time,levels = c("naive","IS_6","IS_24"))
meta$sample<- factor(meta$sample,levels = c("Naive_male_rep5","Naive_male_rep4",
                                            "IS_6h_male_rep1","IS_6h_male_rep2",
                                            "IS_24h_male_rep1","IS_24h_male_rep2"))
#堆积图可视化
library(ggplot2)
library(dplyr)
library(scales)
library(ggsignif)
ggplot(meta, aes(x = time, fill = annotation1)) +
  geom_bar(position = "fill") +   # 每个柱子 0–1，自动变比例
  geom_text(aes(label = ifelse(after_stat(count / sum(count)) > 0.01,  # 可选：比例<1%时不显示，避免拥挤
                               percent(after_stat(count / sum(count)), 
                                       accuracy = 0.1,  # 控制有效数字显示
                                       scale = 100, 
                                       suffix = "%",
                                       signif = 3),   # 关键：保留3位有效数字
                               "")),
            stat = "count",
            position = position_fill(vjust = 0.5),
            colour = "white",
            size = 2.5,
            fontface = "bold") +
  scale_y_continuous(labels = scales::percent_format()) +
  ylab("Cell proportion") +
  xlab("time") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
# annotate("text", x = Inf, y = Inf, label = "* Myelinating Schwann cells (p < 0.05)", 
#hjust = 1.1, vjust = 1.5, size = 4, color = "red")+


#计算各类细胞百分比
library(Seurat)
library(dplyr)
library(ggplot2)
library(tidyverse)
df <- data.frame(ident = seurat$orig.ident,celltype = seurat$annotation1)
count_table <- table(df$ident, df$celltype)
prop_table <- prop.table(count_table, margin = 1) * 100
print(round(prop_table, 2))



#对所有细胞类型统计检验
perform_test <- function(count_table, cell_type, prop_table) {
  idx <- which(colnames(count_table) == cell_type)
  
  # 构建2x2列联表
  test_mat <- matrix(
    c(count_table[1, idx], sum(count_table[1, -idx]),
      count_table[2, idx], sum(count_table[2, -idx])),
    nrow = 2, byrow = TRUE
  )
  
  fisher_res <- fisher.test(test_mat)
  chi_res <- chisq.test(test_mat)
  
  data.frame(
    CellType = cell_type,
    IS_count = count_table["IS", cell_type],
    Naive_count = count_table["Naive", cell_type],
    IS_prop = prop_table["IS", cell_type],
    Naive_prop = prop_table["Naive", cell_type],
    OddsRatio = fisher_res$estimate,
    Fisher_pvalue = fisher_res$p.value,
    Chi2_pvalue = chi_res$p.value
  )
}
stat_results <- do.call(rbind, lapply(colnames(count_table), function(ct) {
  perform_test(count_table, ct, prop_table)
}))
stat_results$p_adj <- p.adjust(stat_results$Fisher_pvalue, method = "BH")#FDR校正

# 添加显著性标记和log2FC
stat_results <- stat_results %>%
  mutate(
    Significance = case_when( p_adj < 0.001 ~ "***",p_adj < 0.01 ~ "**",p_adj < 0.05 ~ "*",TRUE ~ "ns"),
    log2FC = log2((IS_prop + 0.01) / (Naive_prop + 0.01))) %>%
  arrange(p_adj)

plot_df <- data.frame(CellType = rep(colnames(count_table), each = nrow(count_table)),
                      Group = rep(rownames(count_table), times = ncol(count_table)),
                      Count = as.vector(count_table))
plot_df <- plot_df %>%
  group_by(Group) %>%
  mutate(Total = sum(Count),Proportion = Count / Total * 100) %>%
  ungroup() %>%
  left_join(stat_results %>% 
              dplyr::select(CellType, p_adj, Significance, log2FC), by = "CellType")
sig_labels <- plot_df %>%
  group_by(CellType) %>%
  summarise(
    max_prop = max(Proportion),
    Significance = Significance[1],
    .groups = "drop"
  )

p1 <- ggplot(plot_df, aes(x = CellType, y = Proportion, fill = Group)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), 
           width = 0.7, color = "black", linewidth = 0.3) +
  geom_text(data = sig_labels,
            aes(x = CellType, y = max_prop + 2, label = Significance),
            inherit.aes = FALSE, size = 4, fontface = "bold") +
  scale_fill_manual(values = c("Naive" = "#4DBBD5", "IS" = "#E64B35")) +
  labs(title = "Cell Type Proportion Comparison (IS vs Naive)",
       x = "Cell Type", y = "Proportion (%)", fill = "Group") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        plot.title = element_text(hjust = 0.5, face = "bold"),
        legend.position = "top")
ggsave(p1,filename="cell.proportion_comparison(IS.Naive).pdf",width = 12,height = 8)











#差异分析####
#Drd基因集打分
rm(list = ls())
gc()
library(Seurat)
library(AUCell)
library(msigdbr)
library(dplyr)
library(ggplot2)
library(gplots)
library(GSVA)
seurat <- readRDS("/s1/wuqing/PBS_IS/outputs/IS_naive.6.24/5.seurat_annotation1.rds")
table(seurat$annotation1,seurat$time)
DimPlot(seurat,group.by = "annotation1",label = T)
seurat_satglia.c1 <- subset(seurat, subset = annotation1 %in% c("satglia_c1"))
seurat_vascular <- subset(seurat, subset = annotation1 %in% c("vascular"))
seurat_fibro <- subset(seurat, subset = annotation1 %in% c("fibroblast_c1","fibroblast_c2","fibroblast_c3"))
DimPlot(seurat_satglia.c1,label = T,group.by = "time")
table(seurat_vascular$time,seurat_vascular$annotation1)
Idents(seurat_vascular) <-  seurat_vascular$orig.ident
Idents(seurat_satglia.c1)<-  seurat_satglia.c1$orig.ident
Idents(seurat_fibro)<-  seurat_fibro$orig.ident
Idents(seurat)<-  seurat$orig.ident
allMarkers <- FindMarkers(seurat_fibro,
                          only.pos = F,          
                          logfc.threshold = 0.25,
                          min.pct = 0.1,
                          ident.1 = "IS",
                          ident.2 = "Naive") 
allMarkers_filter <- allMarkers %>%
  filter(p_val_adj < 0.05) 
allMarkers_filter <- allMarkers %>%
  filter(pct.1 > 0.1 & p_val_adj < 0.05) %>%
  filter(abs(avg_log2FC)>0.5)
#火山图
library(ggrepel)
library(dplyr)
library(tibble)
colnames(allMarkers_filter)
table(abs(allMarkers_filter$avg_log2FC)>2)
plotdt <- allMarkers_filter %>%
  tibble::rownames_to_column(var = "gene") %>%
  mutate(gene=ifelse(abs(avg_log2FC)>=2,gene,NA))
ggplot(plotdt,aes(x=avg_log2FC,y=-log10(p_val_adj),
                  size=pct.1,
                  color=avg_log2FC))+
  geom_point()+
  ggtitle(label="IS_Naive",subtitle="DEGs")+
  geom_text_repel(aes(label=gene),size=3,color="black")+
  theme_bw()+
  theme(plot.title = element_text(face = "bold",hjust = 0.5),
        plot.background = element_rect(fill = "transparent",color=NA))+
  scale_edge_colour_gradient2(low = "olivedrab",high = "salmon2",
                              mid = "grey",midpoint = 0)+
  scale_size(range = c(1,3))

allMarkers_top20 <- allMarkers %>%
  group_by(cluster) %>%
  slice_max(n = 20, order_by = avg_log2FC)  # 按 log2FC 排序取 top20
allMarkers_top20$gene[allMarkers_top20$cluster==5]

#KEGG/GO####
library(clusterProfiler)
library(org.Mm.eg.db)
library(enrichplot)
allMarkers_filter$gene <- rownames(allMarkers_filter)
ids=bitr(allMarkers_filter$gene,"SYMBOL","ENTREZID","org.Mm.eg.db")
allMarkers_filter=merge(allMarkers_filter,ids,by.x="gene",by.y="SYMBOL")
head(allMarkers_filter)
#将基因按照avg_log2FC的大小降序排列
allMarkers_filter <- allMarkers_filter[order(allMarkers_filter$avg_log2FC,decreasing = T),]
allMarkers_filter_list <- as.numeric(allMarkers_filter$avg_log2FC)
names(allMarkers_filter_list) <- allMarkers_filter$ENTREZID
head(allMarkers_filter_list)

IS.Naive <- names(allMarkers_filter_list)[abs(allMarkers_filter_list)>1]
head(IS.Naive)
#GO
IS.Naive_ego <- enrichGO(IS.Naive,OrgDb = "org.Mm.eg.db",ont = "all",readable = T)
head(IS.Naive_ego)
dotplot(IS.Naive_ego,showCategory=15,title="IS.Naive_GO")
#KEGG
IS.Naive_ekg <- enrichKEGG(IS.Naive,pvalueCutoff = 0.05,organism = "mmu")
head(IS.Naive_ekg)
dotplot(IS.Naive_ekg,showCategory=14,title="IS.Naive_KEGG")








####GO可视化####
GO_result <- read.table("/s1/wuqing/PBS_IS/outputs/IS_naive.6.24/01_组间差异分析_IS/neuron_subtype_DEGs/Go_result_NP.Mrgprd.down.csv",
                        header = T,row.names = 1,sep = ",")
pathways_keep <- c("regulation of synapse structure or activity",
                   "positive regulation of dephosphorylation",
                   "regulation of metal ion transport",
                   "beta-catenin destruction complex" ,
                   "NADH dehydrogenase complex",
                   "voltage-gated potassium channel complex")
plot_data_filtered <- GO_result %>%
  filter(Description %in% pathways_keep) %>%
  mutate(log10_padj = -log10(p.adjust))
range(plot_data_filtered$log10_padj)
ggplot(plot_data_filtered, aes(x = log10_padj, y = reorder(Description, log10_padj))) +
  geom_bar(aes(fill = p.adjust), stat = "identity", width = 0.7) +
  scale_fill_gradient(low = "red", high = "blue", name = "P.adjust") +
  labs(x = "-log10(P.adj)",y=NULL,title = "GO Enrichment of NP_Mrgprd+ downDEGs") +
  theme_bw() +
  theme(axis.text.y = element_text(size = 11, color = "black"),
        axis.text.x = element_text(size = 10, color = "black"),
        axis.title = element_text(size = 12, face = "bold"),
        plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
        legend.position = "right")+
  scale_x_continuous(limits = c(0, 3), breaks = seq(0,3,0.5),expand = c(0, 0))     




GO_result <- read.table("/s1/wuqing/PBS_IS/outputs/IS_naive.6.24/01_组间差异分析_IS/neuron_subtype_DEGs/Go_result_NP.Mrgprd.up.csv",
                        header = T,row.names = 1,sep = ",")
pathways_keep <- c("regulation of synapse structure or activity",
                   "positive regulation of dephosphorylation",
                   "regulation of metal ion transport",
                   "beta-catenin destruction complex" ,
                   "NADH dehydrogenase complex",
                   "voltage-gated potassium channel complex")
plot_data_filtered <- GO_result %>%
  filter(Description %in% pathways_keep) %>%
  mutate(log10_padj = -log10(p.adjust))
range(plot_data_filtered$log10_padj)
ggplot(plot_data_filtered, aes(x = log10_padj, y = reorder(Description, log10_padj))) +
  geom_bar(aes(fill = p.adjust), stat = "identity", width = 0.7) +
  scale_fill_gradient(low = "red", high = "blue", name = "P.adjust") +
  labs(x = "-log10(P.adj)",y=NULL,title = "GO Enrichment of NP_Mrgprd+ downDEGs") +
  theme_bw() +
  theme(axis.text.y = element_text(size = 11, color = "black"),
        axis.text.x = element_text(size = 10, color = "black"),
        axis.title = element_text(size = 12, face = "bold"),
        plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
        legend.position = "right")+
  scale_x_continuous(limits = c(0, 3), breaks = seq(0,3,0.5),expand = c(0, 0))     





#pseudotime-fibro####
rm(list = ls())#清理环境
gc()#释放内存
options(bitmapType='cairo')#图片显示
.libPaths("/home/wuqing/anaconda3/envs/R_4.4/lib/R/library")
library(Seurat)
library(dplyr)
library(tidyverse)
library(Matrix) 
library(patchwork)
library(data.table) 
library(ggplot2)
library(ggrepel)
library(pheatmap)
library(stringr)
library(tibble)
library(monocle)
library(ROGUE)
seurat <- readRDS("/s1/wuqing/PBS_IS/outputs/IS_naive.6.24/5.seurat_annotation1.rds")
seurat <- subset(seurat, subset = annotation1 %in% c("fibroblast_c1","fibroblast_c2","fibroblast_c3"))
#1.构建monocle对象
data <- GetAssayData(seurat,assay="RNA",slot="counts")
colnames(seurat@meta.data)
pd <- new("AnnotatedDataFrame",data=seurat@meta.data[,c("orig.ident","annotation1",
                                                        "seurat_clusters")])
fData <- data.frame(gene_short_name=rownames(data),row.names = row.names(data))#获取基因名
fd <- new("AnnotatedDataFrame",data=fData)
mycds <- newCellDataSet(as(as.matrix(data),"sparseMatrix"),
                        phenoData = pd,
                        featureData = fd,
                        lowerDetectionLimit = 0.5,
                        expressionFamily = negbinomial.size())

#2.数据预处理
mycds <- estimateSizeFactors(mycds)
mycds <- estimateDispersions(mycds,cores=8)
#3.根据基因表达量，离散程度等选取排序基因，2000左右.类似FindVarible
disp_table <- dispersionTable(mycds)
order.genes <- subset(disp_table,mean_expression>=0.005&dispersion_empirical>=
                        1 * dispersion_fit)%>%pull(gene_id)%>%as.character()
#4.排序高变基因
mycds <- setOrderingFilter(mycds,order.genes)
plot_ordering_genes(mycds)
#p <- plot_ordering_genes(mycds)
#ggsave("0.orderGenes.pdf",p,width = 8,height = 6)

#5.初步降维排序
mycds <- reduceDimension(mycds,max_components = 2,reduction_method = "DDRTree")
#  residualModelFormulaStr = "~sample")mycds <- orderCells(mycds)
#6.拟时序排序
mycds <- orderCells(mycds)#root_state=5人工设置起始点


#7.结果可视化
plot_cell_trajectory(mycds,color_by="State")
plot_cell_trajectory(mycds,color_by="Pseudotime")
plot_cell_trajectory(mycds,color_by="annotation1")
plot_cell_trajectory(mycds,color_by="orig.ident")
#树状图
plot_complex_cell_trajectory(mycds,x=1,y=2,
                             color_by = "annotation1")
#细胞密度图
ggplot(pData(mycds),aes(Pseudotime,colour=annotation,fill=celltype))+
  geom_density(bw=0.5,size=1,alpha=0.5)+
  theme_classic()
#指定基因在不同时间表达变化情况
genes <- c(order.genes)[1:4]
plot_genes_in_pseudotime(mycds[genes],color_by="orig.ident")
plot_genes_in_pseudotime(mycds[genes],color_by="annotation")
plot_genes_in_pseudotime(mycds[genes],color_by="Pseudotime")
plot_genes_jitter(mycds[genes],grouping = "State",color_by="annotation")
plot_genes_violin(mycds[genes],grouping = "State",color_by="annotation")
plot_cell_trajectory(mycds,color_by="G3BP2")+
  scale_color_continuous(type = "viridis")

#保存至seurat对象
pdata <- Biobase::pData(mycds)
seurat <- AddMetaData(seurat,metadata=pdata[,c("Pseudotime","State")])
saveRDS(seurat,file="/s1/wuqing/PBS_IS/outputs/IS_naive.6.24/6.seurat_fibro_pseudotime.rds")
#寻找拟时序差异基因
library(BiocGenerics)
library(monocle)
Time_diff <- differentialGeneTest(mycds,
                                  fullModelFormulaStr="~sm.ns(Pseudotime)")
Time_genes <- Time_diff[order(Time_diff$qval),"gene_short_name"][1:200]
write.csv(Time_diff,"time_diff_genes.csv")
plot_pseudotime_heatmap(mycds[Time_genes,],num_clusters=3,show_rownames=T,return_heatmap=T)
saveRDS(mycds,file="/s1/wuqing/PBS_IS/outputs/IS_naive.6.24/7.seurat_fibro_pseudotime_DEGs.rds")



#time_diff分模块
mycds <- read_rds("/s1/wuqing/PBS_IS/outputs/IS_naive.6.24/7.seurat_fibro_pseudotime_DEGs.rds")
modules <- read.csv("module_genes.csv")
Time_diff <- read.csv("time_diff_genes.csv")
Time_diff_filter <- Time_diff[Time_diff$qval<0.05,]%>%arrange(qval)
# 将细胞按pseudotime分为早期和晚期
pseudotime <- mycds$Pseudotime
expr_matrix <- GetAssayData(seurat, slot = "data")
early_cells <- names(pseudotime)[pseudotime < quantile(pseudotime, 0.25)]
late_cells <- names(pseudotime)[pseudotime > quantile(pseudotime, 0.75)]
# 计算log2FC
calc_fc <- function(gene) {
  early_mean <- mean(expr_matrix[gene, early_cells])
  late_mean <- mean(expr_matrix[gene, late_cells])
  log2((late_mean + 0.01) / (early_mean + 0.01))}
Time_diff_filter$log2FC <- sapply(Time_diff_filter$gene_short_name, calc_fc)
# 筛选末期上调基因
late_up_genes <- Time_diff_filter %>%
  filter(log2FC > 0.5, qval < 0.05) %>%
  arrange(desc(log2FC))




#pseudotime-neuron####
rm(list = ls())
gc()
options(bitmapType='cairo')
.libPaths("/home/wuqing/anaconda3/envs/R_4.4/lib/R/library")
library(Seurat)
library(dplyr)
library(tidyverse)
library(Matrix) 
library(patchwork)
library(data.table) 
library(ggplot2)
library(ggrepel)
library(pheatmap)
library(stringr)
library(tibble)
library(monocle)
library(ROGUE)
seurat <- readRDS("/s1/wuqing/PBS_IS/outputs/IS_naive.6.24/6.seurat_neuron_pseudotime.rds")
#1.构建monocle对象
data <- GetAssayData(seurat,assay="RNA",slot="counts")
colnames(seurat@meta.data)
pd <- new("AnnotatedDataFrame",data=seurat@meta.data[,c("orig.ident","subtype","seurat_clusters")])
fData <- data.frame(gene_short_name=rownames(data),
                    row.names = row.names(data))
fd <- new("AnnotatedDataFrame",data=fData)
mycds <- newCellDataSet(as(as.matrix(data),"sparseMatrix"),
                        phenoData = pd,
                        featureData = fd,
                        lowerDetectionLimit = 0.5,
                        expressionFamily = negbinomial.size())
#2.数据预处理
mycds <- estimateSizeFactors(mycds)
mycds <- estimateDispersions(mycds,cores=8)
#3.根据基因表达量，离散程度等选取排序基因，2000左右.类似FindVariblegenes
disp_table <- dispersionTable(mycds)
order.genes <- subset(disp_table,mean_expression>=0.005&dispersion_empirical>=
                        1 * dispersion_fit)%>%pull(gene_id)%>%as.character()
#4.排序高变基因
mycds <- setOrderingFilter(mycds,order.genes)
plot_ordering_genes(mycds)#可视化
#5.初步降维排序
mycds <- reduceDimension(mycds,max_components = 2,reduction_method = "DDRTree")
#6.拟时序排序
mycds <- orderCells(mycds)#root_state=5人工设置起始点
class(mycds)
saveRDS(mycds,"./pseudotime/pseudotime_result.rds")
#7.结果可视化
plot_cell_trajectory(mycds,color_by="State")
plot_cell_trajectory(mycds,color_by="Pseudotime")
plot_cell_trajectory(mycds,color_by="seurat_clusters")
plot_cell_trajectory(mycds,color_by="subtype")
plot_cell_trajectory(mycds,color_by="orig.ident")
plot_genes_in_pseudotime(
        mycds[c("Fbxw7"), ],
        color_by = "subtype",     # 按组别着色
        min_expr = 0.1)
plot_complex_cell_trajectory(mycds,x=1,y=2,color_by = "subtype")#树状图
#细胞密度图
ggplot(pData(mycds),aes(Pseudotime,colour=subtype))+
  geom_density(bw=0.5,size=1,alpha=0.5)+
  theme_classic()
#指定基因在不同时间表达变化情况
genes <- c(order.genes)[1:4]
plot_genes_in_pseudotime(mycds["Klf7",],color_by="subtype")
plot_genes_in_pseudotime(mycds["Klf7",],color_by="Pseudotime")
plot_genes_jitter(mycds[genes],grouping = "State",color_by="subtype")
plot_genes_violin(mycds["Klf7"],grouping = "State",color_by="subtype")
pdata <- Biobase::pData(mycds)
#寻找拟时序差异基因
library(BiocGenerics)
library(monocle)
Time_diff <- read.table("/s1/wuqing/PBS_IS/outputs/IS_naive.6.24/pseudotime/time_diff_genes_neuron.csv",
                      row.names = 1,header = T,sep = ",")
Time_genes <- Time_diff[order(Time_diff$qval),"gene_short_name"][1:200]
plot_pseudotime_heatmap(mycds[Time_genes,],num_clusters=3,show_rownames=T,return_heatmap=T)
#time_diff分模块
Time_diff_filter <- Time_diff[Time_diff$qval<0.05,]%>%arrange(qval)
# 将细胞按pseudotime分为早期和晚期，计算log2FC
pseudotime <- pData(mycds)$Pseudotime
names(pseudotime) <- rownames(pData(mycds))
expr_matrix <- GetAssayData(seurat, slot = "data")
early_cells <- names(pseudotime)[pseudotime < quantile(pseudotime, 0.25)]
late_cells <- names(pseudotime)[pseudotime > quantile(pseudotime, 0.75)]
calc_fc <- function(gene) {
  early_mean <- mean(expr_matrix[gene, early_cells])
  late_mean <- mean(expr_matrix[gene, late_cells])
  log2((late_mean + 0.01) / (early_mean + 0.01))}
Time_diff_filter$log2FC <- sapply(Time_diff_filter$gene_short_name, calc_fc)
# 筛选末期上调基因
late_up_genes <- Time_diff_filter %>%
  filter(log2FC > 0.5, qval < 0.05) %>%
  arrange(desc(log2FC))
write.csv(late_up_genes,"/s1/wuqing/PBS_IS/outputs/IS_naive.6.24/pseudotime/neuron_late_up_genes.csv")
#GSEA富集
library(tidyverse)
library(clusterProfiler)
library(enrichplot)
library(org.Mm.eg.db)
library(dplyr)
library(msigdbr)
late_up_genes <- read.table("/s1/wuqing/PBS_IS/outputs/IS_naive.6.24/pseudotime/neuron_late_up_genes.csv",
                            sep = ",",row.names = 1,header = T)
ids=bitr(late_up_genes$gene_short_name,"SYMBOL","ENTREZID","org.Mm.eg.db")
late_up_genes=merge(late_up_genes,ids,by.x="gene_short_name",by.y="SYMBOL")
#将基因按log2FC降序排列
late_up_genes <- late_up_genes[order(late_up_genes$log2FC,decreasing = T),]
geneList <- as.numeric(late_up_genes$log2FC)
names(geneList) <- as.character(late_up_genes$ENTREZID)
#GSEA
neuron_gseGO <- gseGO(geneList = geneList,
                      OrgDb = org.Mm.eg.db,
                      ont = "All",
                      pvalueCutoff = 0.05)
neuron_gseGO <- neuron_gseGO%>%arrange(desc(NES))
saveRDS(neuron_gseGO,"neuron_gseGO.rds")
neuron_gseGO <- readRDS("./pseudotime/neuron_gseGO.rds")
gseGO_result <- neuron_gseGO@result
gseGO_result <- gseGO_result %>% 
  arrange(desc(NES))
setwd("/s1/wuqing/PBS_IS/outputs/IS_naive.6.24/pseudotime")
write_csv(gseGO_result,"neuron_lateupgenes_gseGO.results_csv")


neuron_gseKEGG <- gseKEGG(geneList = geneList,
                          organism = "mmu",
                          pvalueCutoff = 0.05)
neuron_gseKEGG <- neuron_gseKEGG%>%arrange(desc(NES))
saveRDS(neuron_gseKEGG,"./pseudotime/neuron_gseKEGG.rds")
emapplot(neuron_gseKEGG)
gseKEGG_result <- neuron_gseKEGG@result
gseKEGG_result <- gseKEGG_result %>% 
  arrange(desc(NES))
write_csv(gseKEGG_result,"neuron_lateupgenes_gseKEGG.results_csv")
#可视化
library(enrichplot)
library(ggupset)  
colors <- c("#f7ca64","#43a5bf","#86c697","#a670d6","#ef998a")
gseaplot2(neuron_gseKEGG, 
          geneSetID = c("mmu04080","mmu04750","mmu04270","mmu04082"),
          pvalue_table = TRUE,
          ES_geom = "line")
gseaplot2(neuron_gseGO, 
          geneSetID = c("GO:0098660","GO:0005216","GO:0004888","GO:0004930","GO:0007204",
                        "GO:0001653","GO:0099171","GO:0014059","GO:0048265","GO:0050951",
                        "GO:0048266","GO:0045814"),
          title = "GSEA_enrichment",
          pvalue_table = TRUE,
          ES_geom = "line")
gseaplot2(neuron_gseGO, 
          geneSetID = 25:26,
          pvalue_table = TRUE,
          rel_heights = c(1.5, 0.5, 1),
          subplots = 1:6
          )


