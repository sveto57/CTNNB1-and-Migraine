#cell proportion####
rm(list = ls())
gc()
setwd("/s1/wuqing/PBS_IS/outputs/IS_naive.6.24")
seurat <- readRDS("/s1/wuqing/PBS_IS/outputs/IS_naive.6.24/5.neuron_subtype.rds")
class(seurat$subtype)
table(seurat$subtype)
table(seurat$orig.ident,seurat$subtype)
table(seurat$time,seurat$subtype)
meta <- FetchData(seurat,vars = c("orig.ident", "subtype","time","sample"))


#subtype_marker
table(seurat$annotation1)
table(seurat$subtype)
seurat$subtype <- factor(seurat$subtype, 
                         levels = c("NP_Mrgprd-","NP_Mrgprd+","PEP", "NF_PEP","NF","cLTMR","SST"))
Idents(seurat) <- seurat$subtype

markers <- c("Cd55","Scn11a","Mrgprd",#NP
             "Gal","Tac1",#PEP
             "Nefh","Hapln4",#NF
             "Fam19a4",#cLTMR
             "Sst"#SST
             )
DotPlot(seurat,features = markers,group.by = "subtype")
DEG <- FindAllMarkers(seurat,min.pct=0.25)












meta$subtype <- factor(meta$subtype,
                       levels = c("PEP","NP_Mrgprd-", "NP_Mrgprd+",
                                  "NF_PEP",
                                  "cLTMR","NF",
                                  "SST"))
                                   #   "fibroblast_c1","fibroblast_c2","fibroblast_c3","vascular",
                                    #  "oligodendrocyte", "satglia_c1", "satglia_c2","schwann"
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
ggplot(meta, aes(x = time, fill = subtype)) +
  geom_bar(position = "fill") +   # 每个柱子 0–1，自动变比例
  
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
df <- data.frame(ident = seurat$time,celltype = seurat$subtype)#设置变量
count_table <- table(df$ident, df$celltype)
prop_table <- prop.table(count_table, margin = 1) * 100
print(round(prop_table, 2))
#整体卡方检验（判断整体构成比是否有差异）
library(ggplot2)
library(reshape2)
chisq_test <- chisq.test(count_table)
print(chisq_test)
#针对每个 Subtype 运行独立卡方检验
subtype_tests <- lapply(colnames(prop_table), function(subtype) {
  subtype_count <- count_table[, subtype]
  other_count <- rowSums(count_table) - subtype_count
  contingency <- cbind(subtype_count, other_count)
  test_result <- chisq.test(contingency)
  return(list(subtype = subtype, pvalue = test_result$p.value))})
pval_df <- data.frame(Subtype = sapply(subtype_tests, function(x) x$subtype),
                      pvalue = sapply(subtype_tests, function(x) x$pvalue))
pval_df$label <- ifelse(pval_df$pvalue < 0.001, "***",
                        ifelse(pval_df$pvalue < 0.01, "**",
                               ifelse(pval_df$pvalue < 0.05, "*","ns")))
#可视化数据准备
plot_data <- melt(prop_table,varnames = c("Time", "Subtype"), value.name = "Percentage")
plot_data$Time <- factor(plot_data$Time, levels = c("naive", "IS_6", "IS_24"))
# 计算每个亚型的最大高度用于标注位置
max_height <- aggregate(Percentage ~ Subtype, data = plot_data, FUN = max)
pval_df <- merge(pval_df, max_height, by = "Subtype")
pval_df$y_pos <- pval_df$Percentage + max(plot_data$Percentage) * 0.05
p <- ggplot(plot_data, aes(x = Subtype, y = Percentage, fill = Time)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
  geom_text(data = pval_df, aes(x = Subtype, y = y_pos, label = label),
            inherit.aes = FALSE, size = 3, fontface = "italic") +
  labs(title = "Cell Type Composition Across Time",
       x = NULL, y = "Percentage (%)") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        legend.position = "top") +
  scale_fill_manual(values = c("#00A087", "#4DBBD5", "#E64B35"),
                    labels = c("Naive", "IS-6h", "IS-24h")) +
  ylim(0, max(plot_data$Percentage) * 1.15)  # 扩展y轴留出标注空间
p


#添加显著性标记
if(chi_result$p.value < 0.05) {
  p <- p + annotate("text", x = length(unique(plot_data$Subtype))/2, 
                    y = max(plot_data$Percentage) * 1.05, 
                    label = "***", size = 6)}

print(p)

# 6. 如需成对比较（Bonferroni校正）
adjusted_pvalues <- p.adjust(subtype_pvalues, method = "bonferroni")
result_df <- data.frame(
  Subtype = names(subtype_pvalues),
  P_value = subtype_pvalues,
  Adjusted_P = adjusted_pvalues,
  Significant = ifelse(adjusted_pvalues < 0.05, "Yes", "No")
)
print(result_df)










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






















