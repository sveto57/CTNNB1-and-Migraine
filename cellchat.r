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
library(CellChat)
seurat <- readRDS("/s1/wuqing/PBS_IS/outputs/IS_naive.6.24/5.seurat_annotation1.rds")
seurat
DimPlot(seurat,group.by = "sample")
DimPlot(seurat,group.by = "annotation1",split.by = "time")
#seurat数据完成基本流程分析和去批次


#Cellchat流程####
#1.IS——cellchat####
IS_obj <- subset(seurat, subset = orig.ident == "IS")
IS_obj
IS_data <- GetAssayData(IS_obj,assay = "RNA",slot = "data")#取data数据
IS_metadata <- IS_obj@meta.data[,c("annotation1","orig.ident")]#取metadata数据

#创建cellch对象
cellchat <- createCellChat(object=IS_data)
cellchat <- addMeta(cellchat,meta = IS_metadata)
cellchat <- setIdent(cellchat,ident.use = "annotation1")
levels(cellchat@idents)
table(cellchat@idents)
#小鼠受体配体数据库选择
data(CellChatDB.mouse)
showDatabaseCategory(CellChatDB.mouse)
cellchat@DB <- CellChatDB.mouse
dplyr::glimpse(CellChatDB.mouse$interaction)

cellchat <- subsetData(cellchat,features = NULL)
library(future)
future::plan("multisession",workers=10)
cellchat <- identifyOverExpressedGenes(cellchat)#找高表达配受体基因
cellchat <- identifyOverExpressedInteractions(cellchat)#找HVG对应通路
cellchat <- projectData(cellchat,PPI.mouse)#将基因投射到PPI网络
#CellChat分析：计算每信号通路与配受体对互作的通讯概率推断信号通路水平上的通讯概率,时间较长
cellchat <- computeCommunProb(cellchat,raw.use = T)
#raw.use = T使用data数据分析，是默认参数推荐使用
#raw.use = F使用预处理后的PPI数据，可以找到更多配受体相关的信号通路，可能引进非生物学因素
cellchat <- filterCommunication(cellchat,min.cells = 10)#过滤通讯概率小于10个细胞的低频通路
#汇总相关配受体，计算通路水平上的通信概率
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
cellchat <- netAnalysis_computeCentrality(cellchat,slot.name = "netP")#"netP"表示推断的信号通路的细胞通信网络
IS.net <- subsetCommunication(cellchat)#p值默认小于0.05
setwd("/s1/wuqing/PBS_IS/outputs/IS_naive.6.24")
saveRDS(cellchat,"./cellchat/IS_Cellchat.rds")
write_csv(IS.net,"./cellchat/IS.cellchat_raw.useT.csv")













#1.单组可视化####
rm(list = ls())
gc()
setwd("/s1/wuqing/PBS_IS/outputs/IS_naive.6.24")
cellchat <- readRDS("./cellchat/IS_Cellchat.rds")
IS.net <- read_csv("./cellchat/IS.cellchat_raw.useT.csv")
groupSize <- as.numeric(table(cellchat@idents))
#1.互作网络
par(mfrow=c(1,2))
netVisual_circle(cellchat@net$count,
                 vertex.weight = groupSize,
                 weight.scale = T,
                 label.edge = F,
                 title.name = "number of interactions")#互作数量
netVisual_circle(cellchat@net$weight,
                 vertex.weight = groupSize,
                 weight.scale = T,
                 label.edge = F,
                 title.name = "interaction weights/strength")#互作的强度

#2.展示每个亚群作为source的信号传递
dev.off()
mat <- cellchat@net$weight
par(mfrow=c(3,3),mar=c(1,1,1,1))
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0,nrow = nrow(mat),ncol = ncol(mat),dimnames = dimnames(mat))
  mat2[i, ] <- mat[1,]
  netVisual_circle(mat2,
                   vertex.weight = groupSize,
                   weight.scale = T,
                   arrow.width = 0.2,arrow.size = 0.1,
                   title.name = rownames(mat)[i])
}
dev.off()
par(mfrow=c(2,2),mar=c(1,1,1,1))
mat2 <- matrix(0,nrow = nrow(mat),ncol = ncol(mat),dimnames = dimnames(mat))
mat2[c(2,3), ] <- mat[c(2,3),]
netVisual_circle(mat2,
                 vertex.weight = groupSize,
                 weight.scale = T,
                 arrow.width = 0.2,arrow.size = 0.1,
                 edge.weight.max = max(mat),
                 title.name = rownames(mat)[3])
#3.热图
pathway.show <- IS.ner$pathway_name
cellchat@netP$pathways
pathway.show <- "COLLAGEN"#"CADM"
netVisual_heatmap(cellchat,
                  signaling = pathway.show)
#4.气泡图
levels(cellchat@idents)
netVisual_bubble(cellchat,
                 sources.use = c("schwann"),#指定source
                 targets.use = c("fibroblast_c1","fibroblast_c2","fibroblast_c3","neuron","oligodendrocyte",
                                 "satglia_c1","satglia_c2","schwann","vascular"),#指定target
                 remove.isolate = F)

#5.可视化指定通路基因表达
#参与某条信号通路的所有基因在细胞群中的表达情况展示（小提琴图和气泡图）
library(RColorBrewer)
colors <- brewer.pal(8,"Oranges")#颜色自定义
plotGeneExpression(cellchat,signaling = "COLLAGEN")
plotGeneExpression(cellchat,signaling = "COLLAGEN",type = "dot",col=colors)


#6.网络中心性评分
cellchat <- netAnalysis_computeCentrality(cellchat,slot.name = "netP")
netAnalysis_signalingRole_network(cellchat,signaling = pathway.show,
                                  width = 15,height = 6,font.size = 10)
netAnalysis_signalingRole_scatter(cellchat)
netAnalysis_signalingRole_scatter(cellchat,signaling = c("Glutamate"))#指定通路
netAnalysis_signalingRole_heatmap(cellchat,pattern = "outgoing")
netAnalysis_signalingRole_heatmap(cellchat,pattern = "incoming")
netAnalysis_signalingRole_heatmap(cellchat,pattern = "incoming",
                                  signaling = c("SPP1",))












#2.Naive——cellchat####
Naive_obj <- subset(seurat, subset = orig.ident == "Naive")
Naive_obj
Naive_data <- GetAssayData(Naive_obj,assay = "RNA",slot = "data")#取data数据
Naive_metadata <- Naive_obj@meta.data[,c("annotation1","orig.ident")]#取metadata数据
cellchat <- createCellChat(object=Naive_data)
cellchat <- addMeta(cellchat,meta = Naive_metadata)
cellchat <- setIdent(cellchat,ident.use = "annotation1")
levels(cellchat@idents)
table(cellchat@idents)
data(CellChatDB.mouse)
cellchat@DB <- CellChatDB.mouse
dplyr::glimpse(CellChatDB.mouse$interaction)
cellchat <- subsetData(cellchat,features = NULL)
library(future)
future::plan("multisession",workers=10)
cellchat <- identifyOverExpressedGenes(cellchat)#找高表达配受体基因
cellchat <- identifyOverExpressedInteractions(cellchat)#找HVG对应通路
cellchat <- projectData(cellchat,PPI.mouse)#基因投射到PPI网络
cellchat <- computeCommunProb(cellchat,raw.use = T)
cellchat <- filterCommunication(cellchat,min.cells = 10)#过滤通讯概率小于10个细胞的低频通路
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
cellchat <- netAnalysis_computeCentrality(cellchat,slot.name = "netP")#"netP"表示推断的信号通路的细胞通信网络
Naive.net <- subsetCommunication(cellchat)#p值默认小于0.05
setwd("/s1/wuqing/PBS_IS/outputs/IS_naive.6.24")
saveRDS(cellchat,"./cellchat/Naive_Cellchat.rds")
write_csv(Naive.net,"./cellchat/Naive.cellchat_raw.useT.csv")











rm(list = ls())
gc()
setwd("/s1/wuqing/PBS_IS/outputs/IS_naive.6.24")

IS.cellchat <- readRDS("./cellchat/IS_Cellchat_subtype.rds")
Naive.cellchat<- readRDS("./cellchat/Naive_subtype_Cellchat.rds")
IS.net <- read_csv("./cellchat/IS.cellchat_raw.useT.csv")
Naive.net <- read_csv("./cellchat/Naive.cellchat_raw.useT.csv")
object.list <- list(IS=IS.cellchat,Naive=Naive.cellchat)
cellchat <- mergeCellChat(object.list,add.names = names(object.list))
saveRDS(cellchat,"./cellchat/subtype_Cellchat.rds")
#3.多组cellchat分析####
#1.比较两组互作数目
cellchat <- read_rds("./cellchat/subtype_Cellchat.rds")
cellchat@idents
part_types <- c("NP_Mrgprd+", "NF", "PEP", "NF_PEP","cLTMR","SST","NP_Mrgprd-")
cellchat_subset <- subsetCellChat(cellchat, idents.use = part_types)
gg1 <- compareInteractions(cellchat_subset,show.legend = F,group = c(1,2))
gg2 <- compareInteractions(cellchat_subset,show.legend = F,group = c(1,2),measure = "weight")
gg1+gg2
dev.off()
#2.网络图
#红色表示IS组与Naive组对比，互作次数和强度增加，线越粗表示差异越大
par(mfrow=c(1,2),xpd=T)
netVisual_diffInteraction(cellchat,weight.scale = T)
netVisual_diffInteraction(cellchat,weight.scale = T,
                          measure = "weight",comparison = c(1, 2),
                          color.edge = c("#b2182b", "#2166ac"))
#3.热图
par(mfrow=c(1,1))
h1 <- netVisual_heatmap(cellchat)
h2 <- netVisual_heatmap(cellchat,measure = "weight")
h1+h2

#4.气泡图
par(mfrow=c(1,1))
p1 <- netVisual_bubble(cellchat,
                       sources.use = c("NF_PEP"),#指定source
                       #targets.use = c("fibroblast_c1","cLTMR","NP_Mrgprd-","NP_Mrgprd+","PEP","satglia_c1"),
                       comparison = c(1, 2), 
                       remove.isolate = F)
p2 <- netVisual_bubble(cellchat,
                     #  sources.use = c("fibroblast_c1","cLTMR","NP_Mrgprd-","NP_Mrgprd+","PEP","satglia_c1"),
                       targets.use = c("NF_PEP"),#指定targets
                       comparison = c(1, 2), 
                       remove.isolate = F)
                
p1+p2

#5.保守和特异性信号通路的识别和可视化
gg1 <- rankNet(cellchat,mode = "comparison",stacked = T,do.stat = T)
gg2 <- rankNet(cellchat,mode = "comparison",stacked = F,do.stat = T)
gg1+gg2



#差异分析IS.Naive_cellchat####
setwd("/s1/wuqing/PBS_IS/outputs/IS_naive.6.24/cellchat/neuron_subtype_cellchat_plots")
diff.count <- cellchat@net$IS$count-cellchat@net$Naive$count
write.csv(cellchat@net$IS$count,"output_IScounts_subtype.csv",quote = F)
write.csv(cellchat@net$Naive$count,"output_Naivecounts_subtype.csv",quote = F)
library(pheatmap)
pheatmap(diff.count,
         treeheight_row = "0",treeheight_col = "0",#不画树
         cluster_rows = T,
         cluster_cols = T)











