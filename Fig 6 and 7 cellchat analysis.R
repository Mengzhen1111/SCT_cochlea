setwd("E:\\SCT\\data\\投稿\\CB\\补分析\\cellchat")
#软件安装
devtools::install_github("sqjin/CellChat")

##软件导入
library(CellChat)
library(Seurat)
library(tidyverse)
library(NMF)
library(ggalluvial)
library(patchwork)
library(ggplot2)
library(svglite)
library(patchwork)
#数据读入
object1=readRDS("E17_SCT_L_anno.Rds")
object2=readRDS("P8_SCT_Down_anno.Rds")
object3=readRDS("Adult_SCT_Right_anno.Rds")
##创建SCT assay
object1 <- SCTransform(object1, assay = "Spatial", verbose = TRUE, method = "poisson") 
object2 <- SCTransform(object2, assay = "Spatial", verbose = TRUE, method = "poisson") 
object3 <- SCTransform(object3, assay = "Spatial", verbose = TRUE, method = "poisson") 

###空间数据集准备
#矩阵信息
levels(object1)
levels(object2)
levels(object3)
E_sub <- subset(object1,idents=c("Neuron",
                                 "Kolliker's organ","Tympanic border cell","Hair cell",
                                 "Deiter cell/Pillar cell" ))
P8_sub <- subset(object2,idents=c("Neuron",
                                  "Inner border cell/Inner phalangeal cell","Hensen cell","Tympanic border cell",                   
                                  "Hair cell","Deiter cell/Pillar cell" 
))

A_sub <- subset(object3,idents=c("Neuron",
                                 "Hensen cell",
                                 "Tympanic border cell","Hair cell","Deiter cell/Pillar cell"
))


data.input = Seurat::GetAssayData(E_sub, layer = "data", assay = "SCT") 
data.input = Seurat::GetAssayData(P8_sub, layer = "data", assay = "SCT") 
data.input = Seurat::GetAssayData(A_sub, layer = "data", assay = "SCT") 


#meta信息
meta = data.frame(celltype = Idents(E_sub), #名字自定义
                  row.names = names(Idents(E_sub))) # manually create a dataframe consisting of the cell labels
unique(meta$celltype)

meta = data.frame(celltype = Idents(P8_sub), #名字自定义
                  row.names = names(Idents(P8_sub))) # manually create a dataframe consisting of the cell labels
unique(meta$celltype)

meta = data.frame(celltype = Idents(A_sub), #名字自定义
                  row.names = names(Idents(A_sub))) # manually create a dataframe consisting of the cell labels
unique(meta$celltype)

# 空间图像信息

##E17.5
spatial.locs = Seurat::GetTissueCoordinates(E_sub, scale = NULL, 
                                            cols = c("imagerow", "imagecol")) 
# Scale factors and spot diameters 信息 
scale.factors = jsonlite::fromJSON(txt = 
                                     file.path("E:\\SCT\\data\\cellchat\\E17.5", 'project_setting.json'))
scale.factors = list(spot.diameter = 10, spot = scale.factors$spot_diameter_fullres, # these two information are required
                     fiducial = scale.factors$fiducial_diameter_fullres, hires = scale.factors$tissue_hires_scalef, lowres = scale.factors$tissue_lowres_scalef # these three information are not required
)
##P8
spatial.locs = Seurat::GetTissueCoordinates(P8_sub, scale = NULL, 
                                            cols = c("imagerow", "imagecol")) 
# Scale factors and spot diameters 信息 
scale.factors = jsonlite::fromJSON(txt = 
                                     file.path("E:\\SCT\\data\\cellchat\\P8", 'project_setting.json'))
scale.factors = list(spot.diameter = 10, spot = scale.factors$spot_diameter_fullres, # these two information are required
                     fiducial = scale.factors$fiducial_diameter_fullres, hires = scale.factors$tissue_hires_scalef, lowres = scale.factors$tissue_lowres_scalef # these three information are not required
)

##Adult
spatial.locs = Seurat::GetTissueCoordinates(A_sub, scale = NULL, 
                                            cols = c("imagerow", "imagecol")) 
# Scale factors and spot diameters 信息 
scale.factors = jsonlite::fromJSON(txt = 
                                     file.path("E:\\SCT\\data\\cellchat\\Adult", 'project_setting.json'))
scale.factors = list(spot.diameter = 10, spot = scale.factors$spot_diameter_fullres, # these two information are required
                     fiducial = scale.factors$fiducial_diameter_fullres, hires = scale.factors$tissue_hires_scalef, lowres = scale.factors$tissue_lowres_scalef # these three information are not required
)

##创建对象
##E17.5
E_sub$celltype <- Idents(E_sub)
data.input <- GetAssayData(E_sub, assay = "SCT", layer = "data")
meta <- E_sub@meta.data
cellchat <- createCellChat(
  object = data.input,  # 直接传入表达矩阵
  meta = meta,
  group.by = "celltype"
)

##P8
P8_sub$celltype <- Idents(P8_sub)

data.input <- GetAssayData(P8_sub, assay = "SCT", layer = "data")
meta <- P8_sub@meta.data
cellchat <- createCellChat(
  object = data.input,  # 直接传入表达矩阵
  meta = meta,
  group.by = "celltype"
)

##Adult
A_sub$celltype <- Idents(A_sub)

data.input <- GetAssayData(A_sub, assay = "SCT", layer = "data")
meta <- A_sub@meta.data
cellchat <- createCellChat(
  object = data.input,  # 直接传入表达矩阵
  meta = meta,
  group.by = "celltype"
)

##Cellchat 流程
groupSize <- as.numeric(table(cellchat@idents))
cellChatDB <- CellChatDB.mouse
cellChatDB.use <- subsetDB(cellChatDB, search = "Secreted Signaling")
cellchat@DB <- cellChatDB.use
cellchat <- subsetData(cellchat)
future::plan("multicore", workers = 1)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat) #Identify over-expressed ligand-receptor interactions (pairs) within the used CellChatDB
cellchat <- projectData(cellchat, PPI.mouse) 
cellchat <- computeCommunProb(cellchat,raw.use = F,
                              type = "truncatedMean", trim = 0.1) #如果不想用上一步PPI矫正的结果，raw.use = TRUE即可。
# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat <- filterCommunication(cellchat, min.cells = 1)
#推断信号通路水平的细胞通讯网络
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)

df.net <- subsetCommunication(cellchat)
write.csv(df.net, "E_cellchat_Neu_SE_Secreted signaling.csv")
write.csv(df.net, "P8_cellchat_Neu_SE_Secreted signaling.csv")
write.csv(df.net, "A_cellchat_Neu_SE_Secreted signaling.csv")

df.netp <- subsetCommunication(cellchat, slot.name = "netP")
write.csv(df.netp, "E_cellchat_Neu_SE_Secreted signaling_net_pathway.csv")
write.csv(df.netp, "P8_cellchat_Neu_SE_Secreted signaling_net_pathway.csv")
write.csv(df.netp, "A_cellchat_Neu_SE_Secreted signaling_net_pathway.csv")

##保存cellchat对象
saveRDS(cellchat, "E17.5_Neu_SE_cellchat.Rds")
saveRDS(cellchat, "P8_Neu_SE_cellchat.Rds")
saveRDS(cellchat, "Adult_Neu_SE_cellchat.Rds")

##可视化
##E17.5
cellchat <- readRDS("E17.5_Neu_SE_cellchat.Rds")
df_count <- as.data.frame(cellchat@net$count)
write.csv(df_count,"E_cellchat_Neu_SE_interaction_count_matrix.csv",
          row.names = TRUE)

# 导出通讯强度矩阵（每对细胞群之间的总通讯概率）
df_weight <- as.data.frame(cellchat@net$weight)
write.csv(df_weight, 
          file = "E_cellchat_Neu_SE_interaction_weight_matrix.csv",
          row.names = TRUE)

as.numeric(table(cellchat@idents))
groupSize <- as.numeric(table(cellchat@idents))
#可视化聚合的细胞间通信网络
par(mfrow = c(1,1), xpd=TRUE)
levels(cellchat@idents)

##E17.5
##heatmap可视化通讯数量
mat <- cellchat@net$count
pdf("E17.5_neu_SE_cellchat_number_heatmap.pdf", width = 6, height = 5)
Heatmap(mat,col = colorRampPalette(c("navy","white","firebrick3"))(3),
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        show_column_names = T,
        show_row_names = TRUE)
dev.off()

mat <- cellchat@net$weight
pdf("E17.5_neu_SE_cellchat_weight_heatmap.pdf", width = 6, height = 5)
Heatmap(mat,col = colorRampPalette(c("navy","white","firebrick3"))(3),
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        show_column_names = T,
        show_row_names = TRUE)
dev.off()

##Neuron as source
pdf("E17.5_neu_source_net_number_circle.pdf", width = 4, height = 4)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, 
                 label.edge= F, 
                 sources.use = c(1),
                 targets.use = c(2,3,4,5),
                 color.use = c("#A23EA5",
                               "#FFE4B5","#FD8D3C","#CD2027","#576997"),
                 title.name = "Interaction counts/number")
dev.off()
pdf("E17.5_neu_source_net_weight_circle.pdf", width = 4, height = 4)
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, 
                 label.edge= F, 
                 sources.use = c(1),
                 targets.use = c(2,3,4,5),
                 color.use = c("#A23EA5",
                               "#FFE4B5","#FD8D3C","#CD2027","#576997"),
                 title.name = "Interaction weight")
dev.off()

##HC as source
pdf("E17.5_HC_source_net_number_circle.pdf", width = 4, height = 4)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, 
                 label.edge= F, 
                 sources.use = c(4),
                 targets.use = c(1,2,3,5),
                 color.use = c("#A23EA5",
                               "#FFE4B5","#FD8D3C","#CD2027","#576997"),
                 title.name = "Interaction counts/number")
dev.off()
pdf("E17.5_HC_source_net_weight_circle.pdf", width = 4, height = 4)
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, 
                 label.edge= F, 
                 sources.use = c(4),
                 targets.use = c(1,2,3,5),
                 color.use = c("#A23EA5",
                               "#FFE4B5","#FD8D3C","#CD2027","#576997"),
                 title.name = "Interaction counts/number")
dev.off()

##受体配体气泡图
cellchat@netP$pathways
p <- netVisual_bubble(cellchat, 
                 sources.use = c(1,4),targets.use = c(1,2,3,4,5), 
                 sort.by.source = T,
                 remove.isolate = FALSE)
ggsave("E17.5_Neu_HC_out_all_bubble.pdf", p, width = 6, height = 10)
p <- netVisual_bubble(cellchat, 
                      sources.use = c(1,2,3,4,5),targets.use = c(1,4), 
                      sort.by.source = T,
                      remove.isolate = FALSE)
ggsave("E17.5_Neu_HC_in_all_bubble.pdf", p, width = 6, height = 10)

p <- netVisual_bubble(cellchat, sources.use = c(1,4), signaling = "FGF",
                      targets.use = c(1,2,3,4,5), remove.isolate = FALSE
)
ggsave("E17.5_neu_HC_source_bubble_FGF.pdf", p, width = 6, height = 2.5) #髓系对淋巴的调节
p <- netVisual_bubble(cellchat, sources.use = c(1,2,3,4,5), signaling = "FGF",
                      targets.use = c(1,4), remove.isolate = FALSE
)
ggsave("E17.5_neu_HC_target_bubble_FGF.pdf", p, width = 5, height = 2.5) #髓系对淋巴的调节

##P8
cellchat <- readRDS("P8_Neu_SE_cellchat.Rds")
df_count <- as.data.frame(cellchat@net$count)
write.csv(df_count,"P8_cellchat_Neu_SE_interaction_count_matrix.csv",
          row.names = TRUE)

# 导出通讯强度矩阵（每对细胞群之间的总通讯概率）
df_weight <- as.data.frame(cellchat@net$weight)
write.csv(df_weight, 
          file = "P8_cellchat_Neu_SE_interaction_weight_matrix.csv",
          row.names = TRUE)

as.numeric(table(cellchat@idents))
groupSize <- as.numeric(table(cellchat@idents))
#可视化聚合的细胞间通信网络
par(mfrow = c(1,1), xpd=TRUE)
levels(cellchat@idents)

##heatmap可视化通讯数量
mat <- cellchat@net$count
pdf("P8_neu_SE_cellchat_number_heatmap.pdf", width = 6, height = 5)
Heatmap(mat,col = colorRampPalette(c("navy","white","firebrick3"))(3),
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        show_column_names = T,
        show_row_names = TRUE)
dev.off()

mat <- cellchat@net$weight
pdf("P8_neu_SE_cellchat_weight_heatmap.pdf", width = 6, height = 5)
Heatmap(mat,col = colorRampPalette(c("navy","white","firebrick3"))(3),
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        show_column_names = T,
        show_row_names = TRUE)
dev.off()

##Neuron as source
pdf("P8_neu_source_net_number_circle.pdf", width = 4, height = 4)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, 
                 label.edge= F, 
                 sources.use = c(1),
                 targets.use = c(2,3,4,5,6),
                 color.use = c("#A23EA5",
                               "#FFE4B5","#22816A","#FD8D3C","#CD2027","#576997"),
                 title.name = "Interaction counts/number")
dev.off()
pdf("P8_neu_source_net_weight_circle.pdf", width = 4, height = 4)
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, 
                 label.edge= F, 
                 sources.use = c(1),
                 targets.use = c(2,3,4,5,6),
                 color.use = c("#A23EA5",
                               "#FFE4B5","#22816A","#FD8D3C","#CD2027","#576997"),
                 title.name = "Interaction weight")
dev.off()

##HC as source
pdf("P8_HC_source_net_number_circle.pdf", width = 4, height = 4)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, 
                 label.edge= F, 
                 sources.use = c(5),
                 targets.use = c(1,2,3,4,6),
                 color.use = c("#A23EA5",
                               "#FFE4B5","#22816A","#FD8D3C","#CD2027","#576997"),
                 title.name = "Interaction counts/number")
dev.off()
pdf("P8_HC_source_net_weight_circle.pdf", width = 4, height = 4)
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, 
                 label.edge= F, 
                 sources.use = c(5),
                 targets.use = c(1,2,3,4,6),
                 color.use = c("#A23EA5",
                               "#FFE4B5","#22816A","#FD8D3C","#CD2027","#576997"),
                 title.name = "Interaction counts/number")
dev.off()

##受体配体气泡图
cellchat@netP$pathways
p <- netVisual_bubble(cellchat, 
                      sources.use = c(1,5),targets.use = c(1,2,3,4,5,6), 
                      sort.by.source = T,
                      remove.isolate = FALSE)
ggsave("P8_Neu_HC_out_all_bubble.pdf", p, width = 6, height = 12)
p <- netVisual_bubble(cellchat, 
                      sources.use = c(1,2,3,4,5,6),targets.use = c(1,5), 
                      sort.by.source = T,
                      remove.isolate = FALSE)
ggsave("P8_Neu_HC_in_all_bubble.pdf", p, width = 6, height = 12)

p <- netVisual_bubble(cellchat, sources.use = c(1,5), signaling = "FGF",
                      targets.use = c(1,2,3,4,5,6), remove.isolate = FALSE
)
ggsave("P8_neu_HC_source_bubble_FGF.pdf", p, width = 6, height = 5) #髓系对淋巴的调节
p <- netVisual_bubble(cellchat, sources.use = c(1,2,3,4,5,6), signaling = "FGF",
                      targets.use = c(1,5), remove.isolate = FALSE
)
ggsave("P8_neu_HC_target_bubble_FGF.pdf", p, width = 6, height = 5) #髓系对淋巴的调节


###############################分区cellchat##############################
##重命名需要分析的细胞
levels(object1)
ids <- c("Outer structure",
         "Erythrocytes","Macrophages","Endothelial cell","Neutrophils",
         "Fibrocytes","SLg_Fibrocyte","SLb_Fibrocyte",                    
         "Neuron","Schwann cell",
         "Reissner's membrane",             
         "Stria",
         "Outer sulcus cell",
         "cellchat_celltype",
         "cellchat_celltype",             
         "cellchat_celltype",                
         "cellchat_celltype")
names(ids) <- levels(object1)
object1 <- RenameIdents(object1,ids)
object1$cellchat_cluster <- Idents(object1)
pltdf = data.frame(object1@meta.data, x=object1@images[["sample1"]]@coordinates[["col"]],
                   y=900-object1@images[["sample1"]]@coordinates[["row"]])
object1_color <- c("#F7ECFD", "#F5E8FD","#F3E3FD", "#F1DFFC","#EFDAFC",
                   "#87C986","#C5DEBA","#33A02C",#"Fibrocytes","SL_Fibrocytes","LW_fibrocyte"
                   "#A23EA5",#"Neuron",
                   "#0C2C84",#"Schwann cells"
                   "#E95C59",#"Reissner membrane"
                   "#B49F49",#stria
                   "#3A84E6",
                   "#CD2027")
ggplot(pltdf, 
       aes(x, y, color=cellchat_cluster)) + 
  geom_point( size=2) + 
  scale_color_manual(values = object1_color) +
  theme_bw() +
  theme_light(base_size = 15) +
  labs(title = "") +
  theme(plot.title = element_text(hjust = 0.5))

p = ggplot(pltdf, 
           aes(x, y, color=cellchat_cluster)) + 
  geom_point( size=1) + 
  scale_color_manual(values = object1_color) +
  theme_bw() +
  theme_light(base_size = 15) +
  labs(title = "") +
  theme(plot.title = element_text(hjust = 0.5))

##cellchat_cluster

p+geom_rect(aes(xmin=440,xmax=477,ymin=445,ymax=493), #Apex
            fill = NA,linewidth = 1,color="red") +
  geom_rect(aes(xmin=355,xmax=420,ymin=500,ymax=545), #Middle
            fill = NA,linewidth = 1,color="blue") +
  geom_rect(aes(xmin=400,xmax=460,ymin=390,ymax=445), #Basal
            fill = NA,linewidth = 1,color="black")

d1 = pltdf[which(pltdf$x >= 440 & pltdf$x <= 477 & pltdf$y >= 445 & pltdf$y <= 493),]
d2 = pltdf[which(pltdf$x >= 355 & pltdf$x <= 420 & pltdf$y >= 500 & pltdf$y <= 545),]
d3 = pltdf[which(pltdf$x >= 400 & pltdf$x <= 460 & pltdf$y >= 390 & pltdf$y <= 445),]
#for Neuron
d4 = pltdf[which(pltdf$x >= 402 & pltdf$x <= 450 & pltdf$y >= 461 & pltdf$y <= 495),]
d5 = pltdf[which(pltdf$x >= 380 & pltdf$x <= 402 & pltdf$y >= 423 & pltdf$y <= 453),]
d6 = pltdf[which(pltdf$x >= 370 & pltdf$x <= 398 & pltdf$y >= 476 & pltdf$y <= 508),]


cltp = c("cellchat_celltype") ###change
apex = object1[,colnames(object1) %in% rownames(d1)[d1$cellchat_cluster == cltp]]
apex_se_spots <- colnames(apex)

middle = object1[,colnames(object1) %in% rownames(d2)[d2$cellchat_cluster == cltp]]
middle_se_spots <- colnames(middle)

basal = object1[,colnames(object1) %in% rownames(d3)[d3$cellchat_cluster == cltp]]
basal_se_spots <- colnames(basal)


object1$cellchat_region <- object1$cellchat_cluster
##对目标spot进行命名
# 转换为字符型
object1$cellchat_region <- as.character(object1$cellchat_region)
# 赋值
object1$cellchat_region[colnames(object1) %in% apex_se_spots] <- "Apex_cellchat"
object1$cellchat_region[colnames(object1) %in% middle_se_spots] <- "Middle_cellchat"
object1$cellchat_region[colnames(object1) %in% basal_se_spots] <- "Basal_cellchat"


cltp = c("Neuron") ###change
apex = object1[,colnames(object1) %in% rownames(d4)[d4$region_cluster == cltp]]
apex_sgn_spots <- colnames(apex)
middle = object1[,colnames(object1) %in% rownames(d5)[d5$region_cluster == cltp]]
middle_sgn_spots <- colnames(middle)
basal = object1[,colnames(object1) %in% rownames(d6)[d6$region_cluster == cltp]]
basal_sgn_spots <- colnames(basal)

object1$cellchat_region[colnames(object1) %in% apex_sgn_spots] <- "Apex_Neuron"
object1$cellchat_region[colnames(object1) %in% middle_sgn_spots] <- "Middle_Neuron"
object1$cellchat_region[colnames(object1) %in% basal_sgn_spots] <- "Basal_Neuron"
#再转回因子型
object1$cellchat_region <- as.factor(object1$cellchat_region)

saveRDS(object1,"E:\\SCT\\data\\投稿\\CB\\补分析\\E17.5_anno\\E17_SCT_L_anno.Rds")

##P8
levels(object2)
ids <- c("Outer structure",
         "Neutrophils","Chondrocyte","Pericytes","Osteoblasts",
         "SLg_Fibrocyte","SLb_Fibrocyte",                    
         "Neuron","Schwann cell",
         "Reissner's membrane",             
         "Stria","Spiral prominence","Outer sulcus cell","Interdental cell",
         "Inner sulcus cell","cellchat_celltype",
         "cellchat_celltype",
         "cellchat_celltype",
         "cellchat_celltype",                
         "cellchat_celltype")
names(ids) <- levels(object2)
object2 <- RenameIdents(object2,ids)
object2$cellchat_cluster <- Idents(object2)
pltdf = data.frame(object2@meta.data, x=object2@images[["sample1"]]@coordinates[["col"]],
                   y=900-object2@images[["sample1"]]@coordinates[["row"]])

object2_color <- c("#F7ECFD", "#F5E8FD", "#F3E3FD","#F1DFFC","#EFDAFC",
                   "#87C986","#33A02C",
                   "#A23EA5",#"Neuron",
                   "#0C2C84",#"Schwann cells"
                   "#E95C59",#"Reissner membrane"
                   "#B49F49",#stria
                   "#00B9DB",#Root cell,Spindle cells
                   "#DB4C6C","#FFE4B5",
                   "#FD8D3C","#CD2027"
)
ggplot(pltdf, 
       aes(x, y, color=cellchat_cluster)) + 
  geom_point( size=0.6) + 
  scale_color_manual(values = object2_color) +
  theme_bw() +
  theme_light(base_size = 15) +
  labs(title = "") +
  theme(plot.title = element_text(hjust = 0.5))

p = ggplot(pltdf, 
           aes(x, y, color=cellchat_cluster)) + 
  geom_point( size=0.6) + 
  scale_color_manual(values = object2_color) +
  theme_bw() +
  theme_light(base_size = 15) +
  labs(title = "") +
  theme(plot.title = element_text(hjust = 0.5))

##cellchat_celltype

p+geom_rect(aes(xmin=430,xmax=525,ymin=90,ymax=260),  #Apex
            fill = NA,linewidth = 1,color="red") +
  geom_rect(aes(xmin=525,xmax=600,ymin=100,ymax=150), #Middle
            fill = NA,linewidth = 1,color="blue") +
  geom_rect(aes(xmin=525,xmax=600,ymin=250,ymax=305), #Basal
            fill = NA,linewidth = 1,color="black")
##SGN分区
p + geom_rect(aes(xmin=500,xmax=550,ymin=170,ymax=220), #S1
              fill = NA,linewidth = 1,color="red") +
  geom_rect(aes(xmin=560,xmax=610,ymin=140,ymax=180), #S2
            fill = NA,linewidth = 1,color="blue") +
  geom_rect(aes(xmin=580,xmax=620,ymin=210,ymax=265), #S3
            fill = NA,linewidth = 1,color="black")

#数据提取
#for SE
d1 = pltdf[which(pltdf$x >= 430 & pltdf$x <= 525 & pltdf$y >= 90 & pltdf$y <= 260),]
d2 = pltdf[which(pltdf$x >= 525 & pltdf$x <= 600 & pltdf$y >= 100 & pltdf$y <= 150),]
d3 = pltdf[which(pltdf$x >= 525 & pltdf$x <= 600 & pltdf$y >= 250 & pltdf$y <= 305),]
#for Neuron
d4 = pltdf[which(pltdf$x >= 500 & pltdf$x <= 550 & pltdf$y >= 170 & pltdf$y <= 220),]
d5 = pltdf[which(pltdf$x >= 560 & pltdf$x <= 610 & pltdf$y >= 140 & pltdf$y <= 180),]
d6 = pltdf[which(pltdf$x >= 580 & pltdf$x <= 620 & pltdf$y >= 210 & pltdf$y <= 265),]


##提取区域目标细胞的spot
cltp = c("cellchat_celltype") ###change
apex = object2[,colnames(object2) %in% rownames(d1)[d1$cellchat_cluster == cltp]]
apex_se_spots <- colnames(apex)

middle = object2[,colnames(object2) %in% rownames(d2)[d2$cellchat_cluster == cltp]]
middle_se_spots <- colnames(middle)

basal = object2[,colnames(object2) %in% rownames(d3)[d3$cellchat_cluster == cltp]]
basal_se_spots <- colnames(basal)


object2$cellchat_region <- object2$cellchat_cluster
##对目标spot进行命名
# 转换为字符型
object2$cellchat_region <- as.character(object2$cellchat_region)
# 赋值
object2$cellchat_region[colnames(object2) %in% apex_se_spots] <- "Apex_cellchat"
object2$cellchat_region[colnames(object2) %in% middle_se_spots] <- "Middle_cellchat"
object2$cellchat_region[colnames(object2) %in% basal_se_spots] <- "Basal_cellchat"

##提取区域目标细胞的spot
cltp = c("Neuron") ###change
apex = object2[,colnames(object2) %in% rownames(d4)[d4$cellchat_cluster == cltp]]
apex_sgn_spots <- colnames(apex)
middle = object2[,colnames(object2) %in% rownames(d5)[d5$cellchat_cluster == cltp]]
middle_sgn_spots <- colnames(middle)
basal = object2[,colnames(object2) %in% rownames(d6)[d6$cellchat_cluster == cltp]]
basal_sgn_spots <- colnames(basal)

object2$cellchat_region[colnames(object2) %in% apex_sgn_spots] <- "Apex_Neuron"
object2$cellchat_region[colnames(object2) %in% middle_sgn_spots] <- "Middle_Neuron"
object2$cellchat_region[colnames(object2) %in% basal_sgn_spots] <- "Basal_Neuron"
#再转回因子型
object2$cellchat_region <- as.factor(object2$cellchat_region)
saveRDS(object2,"P8_SCT_Down_anno.Rds")

##分别做不同分区的cellchat
##E17.5
Idents(object1) <- object1$cellchat_region
levels(object1)
sublist <- list(Apex=subset(object1,idents=c("Apex_cellchat","Apex_Neuron")),
                Middle=subset(object1,idents=c("Middle_cellchat","Middle_Neuron")),
                Basal=subset(object1,idents=c("Basal_cellchat","Basal_Neuron")))
##创建对象
for (region in names(sublist)){
  data <- sublist[[region]] 
  Idents(data) <- data$celltype
  data.input <- GetAssayData(data, assay = "SCT", layer = "data")
  meta <- data@meta.data
  meta$celltype <- droplevels(meta$celltype)
  cellchat <- createCellChat(
  object = data.input,  # 直接传入表达矩阵
  meta = meta,
  group.by = "celltype")
groupSize <- as.numeric(table(cellchat@idents))
cellChatDB <- CellChatDB.mouse
cellChatDB.use <- subsetDB(cellChatDB, search = "Secreted Signaling")
cellchat@DB <- cellChatDB.use
cellchat <- subsetData(cellchat)
future::plan("multicore", workers = 1)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat) 
cellchat <- projectData(cellchat, PPI.mouse) 
cellchat <- computeCommunProb(cellchat,raw.use = F,
                              type = "truncatedMean", trim = 0.1) #如果不想用上一步PPI矫正的结果，raw.use = TRUE即可。
cellchat <- filterCommunication(cellchat, min.cells = 1)
#推断信号通路水平的细胞通讯网络
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)

saveRDS(cellchat, paste0("E17.5_Neu_SE_cellchat_",region,".Rds"))

}

##P8
Idents(object2) <- object2$cellchat_region
levels(object2)
sublist <- list(Apex=subset(object2,idents=c("Apex_cellchat","Apex_Neuron")),
                Middle=subset(object2,idents=c("Middle_cellchat","Middle_Neuron")),
                Basal=subset(object2,idents=c("Basal_cellchat","Basal_Neuron")))
##创建对象
for (region in names(sublist)){
  data <- sublist[[region]] 
  Idents(data) <- data$celltype
  data.input <- GetAssayData(data, assay = "SCT", layer = "data")
  meta <- data@meta.data
  meta$celltype <- droplevels(meta$celltype)
  cellchat <- createCellChat(
    object = data.input,  # 直接传入表达矩阵
    meta = meta,
    group.by = "celltype")
  groupSize <- as.numeric(table(cellchat@idents))
  cellChatDB <- CellChatDB.mouse
  cellChatDB.use <- subsetDB(cellChatDB, search = "Secreted Signaling")
  cellchat@DB <- cellChatDB.use
  cellchat <- subsetData(cellchat)
  future::plan("multicore", workers = 1)
  cellchat <- identifyOverExpressedGenes(cellchat)
  cellchat <- identifyOverExpressedInteractions(cellchat) 
  cellchat <- projectData(cellchat, PPI.mouse) 
  cellchat <- computeCommunProb(cellchat,raw.use = F,
                                type = "truncatedMean", trim = 0.1) #如果不想用上一步PPI矫正的结果，raw.use = TRUE即可。
  cellchat <- filterCommunication(cellchat, min.cells = 1)
  #推断信号通路水平的细胞通讯网络
  cellchat <- computeCommunProbPathway(cellchat)
  cellchat <- aggregateNet(cellchat)
  
  saveRDS(cellchat, paste0("P8_Neu_SE_cellchat_",region,".Rds"))
  
}

##分区cellchat可视化
##P8
cellchat_list <- c(P8_Apex=readRDS("P8_Neu_SE_cellchat_Apex.Rds"),
                   P8_middle=readRDS("P8_Neu_SE_cellchat_Middle.Rds"),
                   P8_Basal=readRDS("P8_Neu_SE_cellchat_Basal.Rds"))
output_dir <- "cellchat/P8/"
for (region in names(cellchat_list)) {
  cellchat_obj <- cellchat_list[[region]]
  #导出通讯数量矩阵（每对细胞群之间的L-R对数量）
  df_count <- as.data.frame(cellchat_obj@net$count)
  write.csv(df_count, 
            file = paste0(output_dir, region, "_interaction_count_matrix.csv"),
            row.names = TRUE)
  
  # 导出通讯强度矩阵（每对细胞群之间的总通讯概率）
  df_weight <- as.data.frame(cellchat_obj@net$weight)
  write.csv(df_weight, 
            file = paste0(output_dir, region, "_interaction_weight_matrix.csv"),
            row.names = TRUE)
  # 2. 导出配体-受体对水平通讯（net）
  df_net <- subsetCommunication(cellchat_obj)
  write.csv(df_net, 
            file = paste0(output_dir, region, "_LR_communications.csv"), 
            row.names = FALSE)
  
  # 3. 导出信号通路水平通讯（netP）
  df_netp <- subsetCommunication(cellchat_obj, slot.name = "netP")
  write.csv(df_netp, 
            file = paste0(output_dir, region, "_pathway_communications.csv"), 
            row.names = FALSE)
}
cellchat <- mergeCellChat(cellchat_list, add.names = names(cellchat_list), cell.prefix = TRUE)
cellchat@idents <- unlist(cellchat@idents)
cell_types <- levels(cellchat@idents)
##Neuron as source
pdf("P8_neu_source_net_number_circle_Apex.pdf", width = 4, height = 4)
netVisual_circle(cellchat@net$P8_Apex$count, vertex.weight = groupSize, weight.scale = T, 
                 label.edge= F, 
                 sources.use = c(1),
                 targets.use = c(2,3,4,5,6),
                 color.use = c("#A23EA5",
                               "#FFE4B5","#22816A","#FD8D3C","#CD2027","#576997"),
                 title.name = "P8 Apex")
dev.off()
pdf("P8_neu_source_net_number_circle_Middle.pdf", width = 4, height = 4)
netVisual_circle(cellchat@net$P8_middle$count, vertex.weight = groupSize, weight.scale = T, 
                 label.edge= F, 
                 sources.use = c(1),
                 targets.use = c(2,3,4,5,6),
                 color.use = c("#A23EA5",
                               "#FFE4B5","#22816A","#FD8D3C","#CD2027","#576997"),
                 title.name = "P8 Middle")
dev.off()
pdf("P8_neu_source_net_number_circle_Basal.pdf", width = 4, height = 4)
netVisual_circle(cellchat@net$P8_Basal$count, vertex.weight = groupSize, weight.scale = T, 
                 label.edge= F, 
                 sources.use = c(1),
                 targets.use = c(2,3,4,5,6),
                 color.use = c("#A23EA5",
                               "#FFE4B5","#22816A","#FD8D3C","#CD2027","#576997"),
                 title.name = "P8 Basal")
dev.off()

pdf("P8_neu_source_net_weight_circle_Apex.pdf", width = 4, height = 4)
netVisual_circle(cellchat@net$P8_Apex$weight, vertex.weight = groupSize, weight.scale = T, 
                 label.edge= F, 
                 sources.use = c(1),
                 targets.use = c(2,3,4,5,6),
                 color.use = c("#A23EA5",
                               "#FFE4B5","#22816A","#FD8D3C","#CD2027","#576997"),
                 title.name = "Interaction weight")
dev.off()

pdf("P8_neu_source_net_weight_circle_Middle.pdf", width = 4, height = 4)
netVisual_circle(cellchat@net$P8_middle$weight, vertex.weight = groupSize, weight.scale = T, 
                 label.edge= F, 
                 sources.use = c(1),
                 targets.use = c(2,3,4,5,6),
                 color.use = c("#A23EA5",
                               "#FFE4B5","#22816A","#FD8D3C","#CD2027","#576997"),
                 title.name = "Interaction weight")
dev.off()

pdf("P8_neu_source_net_weight_circle_Basal.pdf", width = 4, height = 4)
netVisual_circle(cellchat@net$P8_Basal$weight, vertex.weight = groupSize, weight.scale = T, 
                 label.edge= F, 
                 sources.use = c(1),
                 targets.use = c(2,3,4,5,6),
                 color.use = c("#A23EA5",
                               "#FFE4B5","#22816A","#FD8D3C","#CD2027","#576997"),
                 title.name = "Interaction weight")
dev.off()

#受体配体气泡图
pdf("P8_neu_HC_source_bubble_regions.pdf", width = 8, height = 15)

netVisual_bubble(cellchat,comparison = c(1,2,3),
                 sources.use = c(1,5),targets.use = c(1,2,3,4,5,6), 
                 sort.by.source = T,
                 remove.isolate = FALSE)

dev.off()
pdf("P8_neu_HC_target_bubble_regions.pdf", width = 8, height = 15)

netVisual_bubble(cellchat,comparison = c(1,2,3),
                 sources.use = c(1,2,3,4,5,6),targets.use = c(1,5), 
                 sort.by.source = T,
                 remove.isolate = FALSE)

dev.off()
##Neuron sorce & target 
pdf("P8_neu_source_bubble_Apex.pdf", width = 4, height = 6)

netVisual_bubble(cellchat_list[[1]],
                 sources.use = c(1),targets.use = c(2,3,4,5,6), 
                 sort.by.source = T,
                 remove.isolate = FALSE)

dev.off()

pdf("P8_neu_target_bubble_Apex.pdf", width = 4, height = 3.5)

netVisual_bubble(cellchat_list[[1]],
                 sources.use = c(2,3,4,5,6),targets.use = c(1), 
                 sort.by.source = T,
                 remove.isolate = FALSE)

dev.off()

pdf("P8_neu_source_bubble_Basal.pdf", width = 4, height = 7)

netVisual_bubble(cellchat_list[[3]],
                 sources.use = c(1),targets.use = c(2,3,4,5,6), 
                 sort.by.source = T,
                 remove.isolate = FALSE)

dev.off()

pdf("P8_neu_target_bubble_Basal.pdf", width = 4, height = 6)

netVisual_bubble(cellchat_list[[3]],
                 sources.use = c(2,3,4,5,6),targets.use = c(1), 
                 sort.by.source = T,
                 remove.isolate = FALSE)

dev.off()

##HC sorce & target 
pdf("P8_HC_target_bubble_Apex.pdf", width = 4, height = 7)

netVisual_bubble(cellchat_list[[1]],
                 sources.use = c(1,2,3,4,6),targets.use = c(5), 
                 sort.by.source = T,
                 remove.isolate = FALSE)

dev.off()

pdf("P8_HC_source_bubble_Apex.pdf", width = 4, height = 5)

netVisual_bubble(cellchat_list[[1]],
                 sources.use = c(5),targets.use = c(1,2,3,4,6), 
                 sort.by.source = T,
                 remove.isolate = FALSE)

dev.off()

pdf("P8_HC_target_bubble_Basal.pdf", width = 4, height = 7)

netVisual_bubble(cellchat_list[[3]],
                 sources.use = c(1,2,3,4,6),targets.use = c(5), 
                 sort.by.source = T,
                 remove.isolate = FALSE)

dev.off()

pdf("P8_HC_source_bubble_Basal.pdf", width = 4, height = 6)

netVisual_bubble(cellchat_list[[3]],
                 sources.use = c(5),targets.use = c(1,2,3,4,6), 
                 sort.by.source = T,
                 remove.isolate = FALSE)

dev.off()

##E17.5
cellchat_list <- c(E_Apex=readRDS("E17.5_Neu_SE_cellchat_Apex.Rds"),
                   E_middle=readRDS("E17.5_Neu_SE_cellchat_Middle.Rds"),
                   E_Basal=readRDS("E17.5_Neu_SE_cellchat_Basal.Rds"))
output_dir <- "cellchat/E17.5/"
for (region in names(cellchat_list)) {
  cellchat_obj <- cellchat_list[[region]]
  #导出通讯数量矩阵（每对细胞群之间的L-R对数量）
  df_count <- as.data.frame(cellchat_obj@net$count)
  write.csv(df_count, 
            file = paste0(output_dir, region, "_interaction_count_matrix.csv"),
            row.names = TRUE)
  
  # 导出通讯强度矩阵（每对细胞群之间的总通讯概率）
  df_weight <- as.data.frame(cellchat_obj@net$weight)
  write.csv(df_weight, 
            file = paste0(output_dir, region, "_interaction_weight_matrix.csv"),
            row.names = TRUE)
  # 2. 导出配体-受体对水平通讯（net）
  df_net <- subsetCommunication(cellchat_obj)
  write.csv(df_net, 
            file = paste0(output_dir, region, "_LR_communications.csv"), 
            row.names = FALSE)
  
  # 3. 导出信号通路水平通讯（netP）
  df_netp <- subsetCommunication(cellchat_obj, slot.name = "netP")
  write.csv(df_netp, 
            file = paste0(output_dir, region, "_pathway_communications.csv"), 
            row.names = FALSE)
}

cellchat <- mergeCellChat(cellchat_list, add.names = names(cellchat_list), cell.prefix = TRUE)
cellchat@idents <- unlist(cellchat@idents)
cell_types <- levels(cellchat@idents)
as.numeric(table(cellchat@idents))
groupSize <- as.numeric(table(cellchat@idents))
##Neuron as source
pdf("E17.5_neu_source_net_number_circle_Apex.pdf", width = 4, height = 4)
netVisual_circle(cellchat@net$E_Apex$count, vertex.weight = groupSize, weight.scale = T, 
                 label.edge= F, 
                 sources.use = c(1),
                 targets.use = c(2,3,4,5),
                 color.use = c("#A23EA5",
                               "#FFE4B5","#FD8D3C","#CD2027","#576997"),
                 title.name = "E17.5 Apex")
dev.off()
pdf("E17.5_neu_source_net_number_circle_Middle.pdf", width = 4, height = 4)
netVisual_circle(cellchat@net$E_middle$count, vertex.weight = groupSize, weight.scale = T, 
                 label.edge= F, 
                 sources.use = c(1),
                 targets.use = c(2,3,4,5),
                 color.use = c("#A23EA5",
                               "#FFE4B5","#FD8D3C","#CD2027","#576997"),
                 title.name = "E17.5 Middle")
dev.off()
pdf("E17.5_neu_source_net_number_circle_Basal.pdf", width = 4, height = 4)
netVisual_circle(cellchat@net$E_Basal$count, vertex.weight = groupSize, weight.scale = T, 
                 label.edge= F, 
                 sources.use = c(1),
                 targets.use = c(2,3,4,5),
                 color.use = c("#A23EA5",
                               "#FFE4B5","#FD8D3C","#CD2027","#576997"),
                 title.name = "E17.5 Basal")
dev.off()

##HC as source
pdf("E17.5_HC_source_net_number_circle_Apex.pdf", width = 4, height = 4)
netVisual_circle(cellchat@net$E_Apex$count, vertex.weight = groupSize, weight.scale = T, 
                 label.edge= F, 
                 sources.use = c(4),
                 targets.use = c(1,2,3,4,5),
                 color.use = c("#A23EA5",
                               "#FFE4B5","#FD8D3C","#CD2027","#576997"),
                 title.name = "E17.5 Apex")
dev.off()
pdf("E17.5_HC_source_net_number_circle_Middle.pdf", width = 4, height = 4)
netVisual_circle(cellchat@net$E_middle$count, vertex.weight = groupSize, weight.scale = T, 
                 label.edge= F, 
                 sources.use = c(4),
                 targets.use = c(1,2,3,4,5),
                 color.use = c("#A23EA5",
                               "#FFE4B5","#FD8D3C","#CD2027","#576997"),
                 title.name = "E17.5 Middle")
dev.off()
pdf("E17.5_HC_source_net_number_circle_Basal.pdf", width = 4, height = 4)
netVisual_circle(cellchat@net$E_Basal$count, vertex.weight = groupSize, weight.scale = T, 
                 label.edge= F, 
                 sources.use = c(4),
                 targets.use = c(1,2,3,4,5),
                 color.use = c("#A23EA5",
                               "#FFE4B5","#FD8D3C","#CD2027","#576997"),
                 title.name = "E17.5 Basal")
dev.off()

p<- rankNet(cellchat, mode = "comparison", stacked = T, comparison = c(1,2,3),do.stat = TRUE,measure = "weight")
ggsave("E17.5_cellchat_region_Compare_pathway_weight.pdf", p, width = 5, height = 6)
p <- rankNet(cellchat, mode = "comparison", stacked = T, comparison = c(1,2,3),do.stat = TRUE,measure = "count")
ggsave("E17.5_cellchat_region_Compare_pathway_count.pdf", p, width = 5, height = 6)

##受体配体气泡图
pdf("E17.5_neu_HC_source_bubble_regions.pdf", width = 6, height = 7)

netVisual_bubble(cellchat,comparison = c(1,2,3),
                 sources.use = c(1,4),targets.use = c(2,3,5), 
                 sort.by.source = T,
                 remove.isolate = FALSE)

dev.off()
pdf("E17.5_neu_HC_target_bubble_regions.pdf", width = 6, height = 7)

netVisual_bubble(cellchat,comparison = c(1,2,3),
                 sources.use = c(2,3,5),targets.use = c(1,4), 
                 sort.by.source = T,
                 remove.isolate = FALSE)

dev.off()

