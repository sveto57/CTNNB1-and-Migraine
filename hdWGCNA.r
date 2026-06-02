#hdWGCNA-fibro####
rm(list = ls())
gc()
library(hdWGCNA)
library(Seurat)
library(tidyverse)
library(WGCNA)
seurat <- readRDS("/s1/wuqing/PBS_IS/outputs/IS_naive.6.24/5.seurat_annotation1.rds")
seurat_fibro <- subset(seurat, subset = annotation1 %in% c("fibroblast_c1","fibroblast_c2","fibroblast_c3"))

DefaultAssay(seurat_fibro) <- "RNA3"
# 1. 准备数据
seurat_fibro <- SetupForWGCNA(seurat_fibro,gene_select = "fraction",fraction = 0.05,wgcna_name = "fibro")
# 2. 构建 metacells
seurat_fibro <- MetacellsByGroups(seurat_fibro,group.by = "annotation1",k = 15,
                                  max_shared = 10,ident.group = "annotation1")
seurat_fibro <- NormalizeMetacells(seurat_fibro)

# 3. 设置表达矩阵
seurat_fibro <- SetDatExpr(seurat_fibro,group.by = "annotation1",
                           group_name = unique(seurat_fibro$annotation1),
                           assay = "RNA3",slot = "data")
# 4. 选择软阈值
seurat_fibro <- TestSoftPowers(seurat_fibro, networkType = "signed")
PlotSoftPowers(seurat_fibro)  # 选择拐点处的值

# 5. 构建网络
seurat_fibro <- ConstructNetwork(seurat_fibro,soft_power = 6,  # 根据上图调整
                                 setDatExpr = FALSE,networkType = "signed",minModuleSize = 30)
PlotDendrogram(seurat_fibro)
# 6. 计算模块特征值和连接度
seurat_fibro <- ModuleEigengenes(seurat_fibro, group.by = "annotation1")
seurat_fibro <- ModuleConnectivity(seurat_fibro)
seurat_fibro <- ResetModuleNames(seurat_fibro, new_name = "M")
# 7. 可视化
seurat_fibro <- ModuleExprScore(seurat_fibro, n_genes = 25)
ModuleFeaturePlot(seurat_fibro, features = "hMEs")
# 8. 提取结果
modules <- GetModules(seurat_fibro)
hub_genes <- GetHubGenes(seurat_fibro, n_hubs = 10)
write.csv(modules, "module_genes.csv")
write.csv(hub_genes, "hub_genes.csv")






# 方法1：模块特征值热图（推荐）
library(ggplot2)
# 获取模块特征值（Module Eigengenes）
MEs <- GetMEs(seurat_fibro, harmonized = TRUE)
seurat_fibro@meta.data <- cbind(seurat_fibro@meta.data, MEs)

# 计算每个亚群的平均模块活性
library(dplyr)
module_avg <- seurat_fibro@meta.data %>%
  group_by(annotation1) %>%
  summarise(across(starts_with("M"), mean)) %>%
  column_to_rownames("annotation1")

module_avg_clean <- module_avg[, colSums(is.na(module_avg)) == 0]
module_avg_clean <- module_avg_clean[rowSums(is.na(module_avg_clean)) == 0, ]
# 绘制热图
library(pheatmap)
pheatmap(
  t(module_avg_clean),
  scale = "row",
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  color = colorRampPalette(c("blue", "white", "red"))(100),
  main = "Module Activity across Cell Types"
)




library(hdWGCNA)
library(clusterProfiler)
library(org.Mm.eg.db)

# 获取模块基因
module2 <- modules %>%filter(module == "M2") %>%pull(gene_name)
gene_ids <- bitr(module2, 
                 fromType = "SYMBOL", 
                 toType = "ENTREZID", 
                 OrgDb = org.Mm.eg.db)

#GO
go_result <- enrichGO(gene = gene_ids$ENTREZID,
                      OrgDb = org.Mm.eg.db,
                      ont = "ALL",
                      pAdjustMethod = "BH",
                      pvalueCutoff = 0.05,
                      qvalueCutoff = 0.2,
                      readable = TRUE )
head(as.data.frame(go_result), 10)
dotplot(go_result, showCategory = 10, split = "ONTOLOGY") +
  facet_grid(ONTOLOGY ~ ., scales = "free")

barplot(go_result, showCategory = 15)
# KEGG
kegg_result <- enrichKEGG(gene = gene_ids$ENTREZID,organism = 'mmu',pvalueCutoff = 0.05)
head(as.data.frame(kegg_result), 10)
dotplot(kegg_result, showCategory = 15)
barplot(kegg_result, showCategory = 15)  

# 可视化
library(ggplot2)
dotplot(go_result, showCategory = 10, split = "ONTOLOGY") 
barplot(kegg_result, showCategory = 15)

















#hdWGCNA-neuron####
rm(list = ls())
gc()
options(bitmapType='cairo')
.libPaths("/home/wuqing/anaconda3/envs/R_4.4/lib/R/library")
library(hdWGCNA)
library(Seurat)
library(tidyverse)
library(WGCNA)
library(patchwork)
setwd("/s1/wuqing/PBS_IS/outputs")
seurat <- readRDS("/s1/wuqing/PBS_IS/outputs/IS_naive.6.24/5.neuron_subtype.rds")
table(seurat$subtype, useNA = "ifany")
DefaultAssay(seurat)
table(seurat$neuron_or_not)
# 1. 设置hdWGCNA对象
seurat <- JoinLayers(seurat)
seurat[["RNA"]] <- as(seurat[["RNA"]], "Assay")
seurat <- SetupForWGCNA(seurat,
                        gene_select = "fraction",#基因选择方式
                        fraction = 0.05,#至少在5%细胞中表达
                        wgcna_name = "hdWGCNA")
# 2. 构建 metacells
seurat <- MetacellsByGroups(seurat,
                            group.by = "subtype",
                            k = 25,
                            max_shared = 10,
                            ident.group = "subtype")
# 减少并行workers
plan("multisession", workers = 8)
options(future.globals.maxSize = 1000 * 1024^2)
seurat <- NormalizeMetacells(seurat)
metacell_obj <- GetMetacellObject(seurat)
Layers(metacell_obj)
metacell_obj <- JoinLayers(metacell_obj)
Layers(metacell_obj)
seurat <- SetMetacellObject(seurat, metacell_obj)
# 3. 设置表达矩阵
seurat <- SetDatExpr(seurat,
                     group.by = "subtype",
                     group_name = unique(seurat$subtype),
                     assay = "RNA",
                     slot ="data")
# 4. 选择软阈值slot参数
seurat <- TestSoftPowers(seurat, networkType = "signed")
plot_list <- PlotSoftPowers(seurat)  # 选择拐点处的值
plot_list[[1]]  # Scale Free R² 图（关键），横线上的第一个数值
plot_list[[2]]  # Mean Connectivity 图
# 5. 构建网络
setwd("/s1/wuqing/PBS_IS/outputs/IS_naive.6.24/hdWGCNA/neuron")
seurat <- ConstructNetwork(seurat,
                           soft_power = 6,  # 根据上图调整
                           setDatExpr = FALSE,
                           networkType = "signed",
                           minModuleSize = 30)
PlotDendrogram(seurat)
# 6. 计算模块特征值和连接度
options(future.globals.maxSize = 4000 * 1024^2)
seurat <- ModuleEigengenes(seurat, group.by = "subtype")
seurat <- ModuleConnectivity(seurat)
seurat <- ResetModuleNames(seurat, new_name = "M")
setwd("/s1/wuqing/PBS_IS/outputs/IS_naive.6.24")
saveRDS(seurat,"9.neuron_subtype_hdWGCNA.rds")
# 7. 可视化
seurat <- ModuleExprScore(seurat, n_genes = 25)
ModuleFeaturePlot(seurat, features = "hMEs")
# 8. 提取结果
modules <- GetModules(seurat)
table(modules$module)
hub_genes <- GetHubGenes(seurat, n_hubs = 10)
table(hub_genes$module)
setwd("/s1/wuqing/PBS_IS/outputs/IS_naive.6.24/hdWGCNA/neuron")
write.csv(modules, "neuron.module_genes.csv")
write.csv(hub_genes, "neuron.hub_genes.csv")





#可视化####
rm(list = ls())
gc()
library(hdWGCNA)
library(Seurat)
library(tidyverse)
library(WGCNA)
library(patchwork)
options(bitmapType='cairo')
.libPaths("/home/wuqing/anaconda3/envs/R_4.4/lib/R/library")
# 热图（推荐）
setwd("/s1/wuqing/PBS_IS/outputs/IS_naive.6.24/hdWGCNA/neuron")
library(ggplot2)
modules <- read.table( "neuron.module_genes.csv",sep = ",",header = T,row.names = 1)
hub_genes <- read.table( "neuron.hub_genes.csv",sep = ",",header = T,row.names = 1)
setwd("/s1/wuqing/PBS_IS/outputs/IS_naive.6.24")
seurat <- readRDS("9.neuron_subtype_hdWGCNA.rds")
MEs <- GetMEs(seurat, harmonized = TRUE)#MEs是模块特征基因
seurat@meta.data <- cbind(seurat@meta.data, MEs)


#1.小提琴图比较各模块在两组的表达
library(reshape2)
MEs$group <- seurat$time
colnames(MEs)
MEs <- MEs[,c("M1","M2","M3","M4","M5","M6","M7","M8","M9","M10","M11","M12","group")]
MEs_long <- melt(MEs, id.vars = "group")
MEs_long$group <- factor(MEs_long$group,levels = c("naive", "IS_6", "IS_24"))
ggplot(MEs_long, aes(x = variable, y = value, fill = group)) +
  geom_violin(trim = FALSE) +
  geom_boxplot(width = 0.1, position = position_dodge(0.9)) +
  theme_classic() +
  labs(x = "Module", y = "Module Eigengene", fill = "Group") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggplot(MEs_long, aes(x = group, y = value, fill = group)) +
  geom_violin(trim = FALSE) +
  geom_boxplot(width = 0.2) +
  facet_wrap(~variable, scales = "free_y") +
  theme_bw() +
  labs(x = "Group", y = "Module Eigengene")

#2.热图可视化模块表达模式
library(ComplexHeatmap)
library(circlize)
#按组计算每个模块的平均表达并标记统计学差异
MEs_avg <- MEs %>%
  group_by(group) %>%
  summarise(across(starts_with("M"), mean)) %>%
  column_to_rownames("group")
MEs_avg <- t(MEs_avg)
MEs_avg <- MEs_avg[,c("naive", "IS_6", "IS_24")]
library(ggpubr)
me_stats <- data.frame()
for(me in colnames(MEs)[colnames(MEs) != "group"]) {
  test_result <- kruskal.test(MEs[[me]] ~ MEs$group)
  me_stats <- rbind(me_stats, 
                    data.frame(module = me, 
                               p_value = test_result$p.value))}
me_stats$p_adj <- p.adjust(me_stats$p_value, method = "BH")
add_significance <- function(module_name, p_adj) {
  if(p_adj < 0.001) return(paste0(module_name, " ***"))
  else if(p_adj < 0.01) return(paste0(module_name, " **"))
  else if(p_adj < 0.05) return(paste0(module_name, " *"))
  else return(module_name)
}
# 为MEs_avg的行名添加显著性标记
rownames_with_sig <- sapply(rownames(MEs_avg), function(me) {
  p_val <- me_stats$p_adj[me_stats$module == me]
  if(length(p_val) > 0) {add_significance(me, p_val)} else {me}})
MEs_avg_annotated <- MEs_avg
rownames(MEs_avg_annotated) <- rownames_with_sig
Heatmap(MEs_avg_annotated,
        name = "ME",
        cluster_rows = TRUE,
        cluster_columns = FALSE,
        col = colorRamp2(c(min(MEs_avg), 0, max(MEs_avg)), 
                         c("blue", "white", "red")),
        row_names_gp = gpar(fontsize = 10),
        column_names_gp = gpar(fontsize = 12))



#3.热图可视化计算各模块在每个细胞亚群的活性
library(dplyr)
duplicated_cols <- duplicated(colnames(seurat@meta.data))
table(duplicated_cols)
meta_unique <- seurat@meta.data[, !duplicated(colnames(seurat@meta.data))]
module_avg <- meta_unique %>%
  group_by(subtype) %>%
  summarise(across(starts_with("M"), mean, .names = "{.col}")) %>%
  column_to_rownames("subtype")
module_avg_clean <- module_avg[, colSums(is.na(module_avg)) == 0]
module_avg_clean <- module_avg_clean[rowSums(is.na(module_avg_clean)) == 0, ]
rownames(module_avg_clean)
colnames(module_avg_clean)
module_avg_clean <- module_avg_clean[c("cLTMR","NF","NF_PEP","PEP","NP_Mrgprd-","NP_Mrgprd+","SST"),
                                     c("M1","M2","M3","M4","M5","M6","M7","M8","M9","M10","M11","M12")]
library(pheatmap)
pheatmap(t(module_avg_clean),
         scale = "row",
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         color = colorRampPalette(c("blue", "white", "red"))(100),
         main = "Module Activity across Cell Types")
  
#4.可视化hub_genes
genes_to_plot <- hub_genes %>% 
  filter(module == "M1") %>% 
  pull(gene_name)
seurat$time <- factor(seurat$time,levels = c("naive","IS_6","IS_24"))
VlnPlot(seurat, features = genes_to_plot[1:6],group.by = "time",ncol = 3)
DotPlot(seurat,features = genes_to_plot,group.by = "time") +
  coord_flip() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))



#6.在降维图上展示模块特征基因的表达
FeaturePlot(seurat, split.by = "time",ncol = 3, 
            features = c("M1", "M4", "M10"))





#富集####
# 获取某个模块基因
library(hdWGCNA)
library(clusterProfiler)
library(org.Mm.eg.db)
module1.4.10<- modules %>%filter(module %in% c("M1","M4","M10")) %>%pull(gene_name)
module1.4.10_hub <- hub_genes %>%filter(module %in% c("M1","M4","M10")) %>%pull(gene_name)
setwd("/s1/wuqing/PBS_IS/outputs/IS_naive.6.24/hdWGCNA/neuron")
write.csv(module1.4.10, "module1.4.10.csv")
write.csv(module1.4.10_hub, "module1.4.10_hub.csv")
library(tidyverse)
library(clusterProfiler)
library(enrichplot)
library(org.Mm.eg.db)
library(dplyr)
library(msigdbr)
gene_ids <- bitr(module1.4.10, 
                 fromType = "SYMBOL", 
                 toType = "ENTREZID", 
                 OrgDb = org.Mm.eg.db)
#GO
go_result <- enrichGO(gene = gene_ids$ENTREZID,
                      OrgDb = org.Mm.eg.db,
                      ont = "ALL",
                      pAdjustMethod = "BH",
                      pvalueCutoff = 0.05,
                      qvalueCutoff = 0.2,
                      readable = TRUE )
Go_result <- as.data.frame(go_result)
setwd("/s1/wuqing/PBS_IS/outputs/IS_naive.6.24/hdWGCNA/neuron")
write.csv(Go_result, "module1.4.10_go.result.csv")

dotplot(go_result, showCategory = 10, split = "ONTOLOGY") +
  facet_grid(ONTOLOGY ~ ., scales = "free")
barplot(go_result, showCategory = 25)
#KEGG
kegg_result <- enrichKEGG(gene = gene_ids$ENTREZID,organism = 'mmu',pvalueCutoff = 0.05)
KEGG_result <- as.data.frame(kegg_result)
write.csv(KEGG_result, "module1.4.10_kegg.result.csv")
barplot(kegg_result, showCategory = 5)  
pathways_keep <- c("Cholinergic synapse - Mus musculus (house mouse)", 
                   "Circadian entrainment - Mus musculus (house mouse)",
                   "Chemokine signaling pathway - Mus musculus (house mouse)",
                   "Dopaminergic synapse - Mus musculus (house mouse)" ,
                   "Aldosterone synthesis and secretion - Mus musculus (house mouse)",#神经类固醇
                   "Adherens junction - Mus musculus (house mouse)",
                   "Cortisol synthesis and secretion - Mus musculus (house mouse)",
                   "Hormone signaling - Mus musculus (house mouse)",
                   "Glutamatergic synapse - Mus musculus (house mouse)",
                   "Endocytosis - Mus musculus (house mouse)",
                   "Rap1 signaling pathway - Mus musculus (house mouse)",
                   "Neuroactive ligand signaling - Mus musculus (house mouse)",
                   "GABAergic synapse - Mus musculus (house mouse)",
                   "Oxytocin signaling pathway - Mus musculus (house mouse)",
                   "Estrogen signaling pathway - Mus musculus (house mouse)",
                   "Long-term depression - Mus musculus (house mouse)",
                   "Growth hormone synthesis, secretion and action - Mus musculus (house mouse)",
                   "MAPK signaling pathway - Mus musculus (house mouse)",
                   "cAMP signaling pathway - Mus musculus (house mouse)",
                   "Phospholipase D signaling pathway - Mus musculus (house mouse)",
                   "Inflammatory mediator regulation of TRP channels - Mus musculus (house mouse)")
plot_data_filtered <- KEGG_result %>%
  filter(Description %in% pathways_keep) %>%
  mutate(log10_padj = -log10(p.adjust))
range(plot_data_filtered$log10_padj)
ggplot(plot_data_filtered, aes(x = log10_padj, y = reorder(Description, log10_padj))) +
  geom_bar(aes(fill = p.adjust), stat = "identity", width = 0.7) +
  scale_fill_gradient(low = "red", high = "blue", name = "P.adjust") +
  labs(x = "-log10(P.adj)",y=NULL,title = "KEGG Enrichment of module1&4&10") +
  theme_bw() +
  theme(axis.text.y = element_text(size = 11, color = "black"),
        axis.text.x = element_text(size = 10, color = "black"),
        axis.title = element_text(size = 12, face = "bold"),
        plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
        legend.position = "right")+
  scale_x_continuous(limits = c(0, 4), breaks = seq(0,4,0.5),expand = c(0, 0))     








