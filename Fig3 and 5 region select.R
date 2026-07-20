library(ggunchull)
library(ggplot2)
library(Seurat)
library(tidyverse)
library(Matrix)
library(ggplot2)
library(ggpubr)
library(gridExtra)
library(reshape2)
library(readr)
library(config)
library(gridExtra)
library(png)
library(patchwork)
library(readxl)
library(patchwork)
library(RColorBrewer)
library(colorRamp2)
library(ComplexHeatmap)
library(semla)
library(ggforce)
library(ggplot2)
library(tidydr)
library(cowplot)
library(scCustomize)
library(Nebulosa)
library(ggVennDiagram)
object1_anno <- readRDS("E17_SCT_R_anno.Rds")
object2_anno <- readRDS("E17_SCT_L_anno.Rds")
object3_anno <- readRDS("P8_SCT_Up_anno.Rds")
object4_anno <- readRDS("P8_SCT_Down_anno.Rds")
object6_anno <- readRDS("Adult_SCT_Right_anno.Rds")


################################   For E17.5 R   ############################
Idents(object1_anno) <- object1_anno$celltype
levels(object1_anno)

ids <- c("Outer structure",
         "Macrophages", "Neutrophils","Erythrocytes","Ube2c+",
         "SLg_Fibrocyte","SLb_Fibrocyte",                    
         "Neuron","Schwann cell",
         "Reissner's membrane",             
         "Stria",                            
         "Outer sulcus cell",
         "Sensory Epithelium",              
         "Sensory Epithelium",                 
         "Sensory Epithelium")
names(ids) <- levels(object1_anno)
object1_anno <- RenameIdents(object1_anno,ids)
object1_anno$region_cluster <- Idents(object1_anno)
pltdf = data.frame(object1_anno@meta.data, x=object1_anno@images[["sample1"]]@coordinates[["col"]],
                   y=900-object1_anno@images[["sample1"]]@coordinates[["row"]])

object1_color <- c("#F7ECFD", "#F5E8FD",
                   "#F3E3FD", "#F1DFFC","#EFDAFC",
                   "#C5DEBA","#33A02C",#"Fibrocytes","SL_Fibrocytes","LW_fibrocyte"
                   "#A23EA5",#"Neuron",
                   "#0C2C84",#"Schwann cells"
                   "#E95C59",#"Reissner membrane"
                   "#B49F49",#stria
                   "#3A84E6",
                   "#CD2027")
ggplot(pltdf, 
           aes(x, y, color=region_cluster)) + 
  geom_point( size=1) + 
  scale_color_manual(values = object1_color) +
  theme_bw() +
  theme_light(base_size = 15) +
  labs(title = "") +
  theme(plot.title = element_text(hjust = 0.5))

p = ggplot(pltdf, 
           aes(x, y, color=region_cluster)) + 
  geom_point( size=1) + 
  scale_color_manual(values = object1_color) +
  theme_bw() +
  theme_light(base_size = 15) +
  labs(title = "") +
  theme(plot.title = element_text(hjust = 0.5))

##Sensory Epithelium

p+geom_rect(aes(xmin=630,xmax=740,ymin=250,ymax=330), #Apex
            fill = NA,linewidth = 1,color="red") +
  geom_rect(aes(xmin=580,xmax=630,ymin=270,ymax=310), #Middle
            fill = NA,linewidth = 1,color="blue") +
  geom_rect(aes(xmin=675,xmax=725,ymin=210,ymax=250), #Basal
            fill = NA,linewidth = 1,color="black")
ggsave("E17.5_SCT_R_Sensory Epithelium_region_select.pdf",width = 10,height = 8)
##SGN分区
p + geom_rect(aes(xmin=625,xmax=660,ymin=285,ymax=325), #Apex
              fill = NA,linewidth = 1,color="red") +
  geom_rect(aes(xmin=665,xmax=705,ymin=260,ymax=295), #Middle
            fill = NA,linewidth = 1,color="blue") +
  geom_rect(aes(xmin=610,xmax=670,ymin=220,ymax=260), #Basal
            fill = NA,linewidth = 1,color="black")
ggsave("E17.5_SCT_R_SGN_region_select.pdf",width = 10,height = 8)

#数据提取
#for SE
d1 = pltdf[which(pltdf$x >= 630 & pltdf$x <= 740 & pltdf$y >= 250 & pltdf$y <= 330),]
d2 = pltdf[which(pltdf$x >= 580 & pltdf$x <= 630 & pltdf$y >= 270 & pltdf$y <= 310),]
d3 = pltdf[which(pltdf$x >= 675 & pltdf$x <= 725 & pltdf$y >= 210 & pltdf$y <= 250),]
#for Neuron
d4 = pltdf[which(pltdf$x >= 625 & pltdf$x <= 660 & pltdf$y >= 285 & pltdf$y <= 325),]
d5 = pltdf[which(pltdf$x >= 665 & pltdf$x <= 705 & pltdf$y >= 260 & pltdf$y <= 295),]
d6 = pltdf[which(pltdf$x >= 610 & pltdf$x <= 670 & pltdf$y >= 220 & pltdf$y <= 260),]


##提取区域目标细胞的spot
cltp = c("Sensory Epithelium") ###change
apex = object1_anno[,colnames(object1_anno) %in% rownames(d1)[d1$region_cluster == cltp]]
apex_se_spots <- colnames(apex)

middle = object1_anno[,colnames(object1_anno) %in% rownames(d2)[d2$region_cluster == cltp]]
middle_se_spots <- colnames(middle)

basal = object1_anno[,colnames(object1_anno) %in% rownames(d3)[d3$region_cluster == cltp]]
basal_se_spots <- colnames(basal)


object1_anno$region_annotation <- object1_anno$region_cluster
##对目标spot进行命名
# 转换为字符型
object1_anno$region_annotation <- as.character(object1_anno$region_annotation)
# 赋值
object1_anno$region_annotation[colnames(object1_anno) %in% apex_se_spots] <- "Apex_SE"
object1_anno$region_annotation[colnames(object1_anno) %in% middle_se_spots] <- "Middle_SE"
object1_anno$region_annotation[colnames(object1_anno) %in% basal_se_spots] <- "Basal_SE"

##提取区域目标细胞的spot
cltp = c("Neuron") ###change
apex = object1_anno[,colnames(object1_anno) %in% rownames(d4)[d4$region_cluster == cltp]]
apex_sgn_spots <- colnames(apex)
middle = object1_anno[,colnames(object1_anno) %in% rownames(d5)[d5$region_cluster == cltp]]
middle_sgn_spots <- colnames(middle)
basal = object1_anno[,colnames(object1_anno) %in% rownames(d6)[d6$region_cluster == cltp]]
basal_sgn_spots <- colnames(basal)

object1_anno$region_annotation[colnames(object1_anno) %in% apex_sgn_spots] <- "Apex_Neuron"
object1_anno$region_annotation[colnames(object1_anno) %in% middle_sgn_spots] <- "Middle_Neuron"
object1_anno$region_annotation[colnames(object1_anno) %in% basal_sgn_spots] <- "Basal_Neuron"
#再转回因子型
object1_anno$region_annotation <- as.factor(object1_anno$region_annotation)

Idents(object1_anno) <- object1_anno$region_annotation
levels(object1_anno)
object1_anno <- subset(object1_anno,idents = c("Outer structure",
                                               "Macrophages", "Neutrophils","Erythrocytes","Ube2c+",
                                               "SLg_Fibrocyte","SLb_Fibrocyte",                    
                                               "Apex_Neuron","Middle_Neuron","Basal_Neuron",
                                               "Schwann cell",
                                               "Reissner's membrane",             
                                               "Stria",                            
                                               "Outer sulcus cell",
                                               "Apex_SE","Middle_SE","Basal_SE"))
Idents(object1_anno) <- factor(Idents(object1_anno), levels =c("Outer structure",
                                                               "Macrophages", "Neutrophils","Erythrocytes","Ube2c+",
                                                               "SLg_Fibrocyte","SLb_Fibrocyte",                    
                                                               "Apex_Neuron","Middle_Neuron","Basal_Neuron",
                                                               "Schwann cell",
                                                               "Reissner's membrane",             
                                                               "Stria",                            
                                                               "Outer sulcus cell",
                                                               "Apex_SE","Middle_SE","Basal_SE"))

object1_anno$region_annotation <- Idents(object1_anno)

pltdf = data.frame(object1_anno@meta.data, x=object1_anno@images[["sample1"]]@coordinates[["col"]],
                   y=900-object1_anno@images[["sample1"]]@coordinates[["row"]])

object1_color <- c("#F5F5F5","#F5F5F5",
                   "#F5F5F5","#F5F5F5","#F5F5F5",
                   "#F5F5F5","#F5F5F5",#"Fibrocytes","SL_Fibrocytes","LW_fibrocyte"
                   "#A23EA5", "#2E86C1","#1ABC9C",#"Neuron",
                   "#F5F5F5",#"Schwann cells"
                   "#F5F5F5",#"Reissner membrane"
                   "#F5F5F5",#stria
                   "#F5F5F5",
                   "#F5F5F5","#F5F5F5","#F5F5F5")
ggplot(pltdf, 
       aes(x, y, color=region_annotation)) + 
  geom_point( size=1) + 
  scale_color_manual(values = object1_color) +
  theme_bw() +
  theme_light(base_size = 15) +
  labs(title = "") +
  theme(plot.title = element_text(hjust = 0.5))

ggsave("E17.5_SCT_R_SGN_region_select_result.pdf",width = 10,height = 8)

saveRDS(object1_anno,"E17_SCT_R_anno.Rds")
remove(object1_anno)

############################### For E17.5 L ###################
Idents(object2_anno) <- object2_anno$celltype
levels(object2_anno)

ids <- c("Outer structure",
         "Erythrocytes","Macrophages","Endothelial cell","Neutrophils",
         "Fibrocytes","SLg_Fibrocyte","SLb_Fibrocyte",                    
         "Neuron","Schwann cell",
         "Reissner's membrane",             
         "Stria",
         "Outer sulcus cell",
         "Sensory Epithelium",
         "Tympanic border cell",             
         "Sensory Epithelium",                
         "Sensory Epithelium")
names(ids) <- levels(object2_anno)
object2_anno <- RenameIdents(object2_anno,ids)
object2_anno$region_cluster <- Idents(object2_anno)
pltdf = data.frame(object2_anno@meta.data, x=object2_anno@images[["sample1"]]@coordinates[["col"]],
                   y=900-object2_anno@images[["sample1"]]@coordinates[["row"]])

object2_color <- c("#F7ECFD", "#F5E8FD","#F3E3FD", "#F1DFFC","#EFDAFC",
                   "#87C986","#C5DEBA","#33A02C",#"Fibrocytes","SL_Fibrocytes","LW_fibrocyte"
                   "#A23EA5",#"Neuron",
                   "#0C2C84",#"Schwann cells"
                   "#E95C59",#"Reissner membrane"
                   "#B49F49",#stria
                   "#3A84E6",
                   "#CD2027",
                   "#F898CB")
ggplot(pltdf, 
       aes(x, y, color=region_cluster)) + 
  geom_point( size=2) + 
  scale_color_manual(values = object2_color) +
  theme_bw() +
  theme_light(base_size = 15) +
  labs(title = "") +
  theme(plot.title = element_text(hjust = 0.5))

p = ggplot(pltdf, 
           aes(x, y, color=region_cluster)) + 
  geom_point( size=1) + 
  scale_color_manual(values = object2_color) +
  theme_bw() +
  theme_light(base_size = 15) +
  labs(title = "") +
  theme(plot.title = element_text(hjust = 0.5))

##Sensory Epithelium

p+geom_rect(aes(xmin=440,xmax=477,ymin=445,ymax=493), #Apex
            fill = NA,linewidth = 1,color="red") +
  geom_rect(aes(xmin=355,xmax=420,ymin=500,ymax=545), #Middle
            fill = NA,linewidth = 1,color="blue") +
  geom_rect(aes(xmin=400,xmax=460,ymin=390,ymax=445), #Basal
            fill = NA,linewidth = 1,color="black")
ggsave("E17.5_SCT_L_Sensory Epithelium_region_select.pdf",width = 10,height = 8)
##SGN分区
p + geom_rect(aes(xmin=402,xmax=450,ymin=461,ymax=495), #S1
              fill = NA,linewidth = 1,color="red") +
  geom_rect(aes(xmin=380,xmax=402,ymin=423,ymax=453), #S2
            fill = NA,linewidth = 1,color="blue") +
  geom_rect(aes(xmin=370,xmax=398,ymin=476,ymax=508), #S3
            fill = NA,linewidth = 1,color="black")
ggsave("E17.5_SCT_L_SGN_region_select.pdf",width = 8.3,height = 6.8)

#数据提取
#for SE
d1 = pltdf[which(pltdf$x >= 440 & pltdf$x <= 477 & pltdf$y >= 445 & pltdf$y <= 493),]
d2 = pltdf[which(pltdf$x >= 355 & pltdf$x <= 420 & pltdf$y >= 500 & pltdf$y <= 545),]
d3 = pltdf[which(pltdf$x >= 400 & pltdf$x <= 460 & pltdf$y >= 390 & pltdf$y <= 445),]
#for Neuron
d4 = pltdf[which(pltdf$x >= 402 & pltdf$x <= 450 & pltdf$y >= 461 & pltdf$y <= 495),]
d5 = pltdf[which(pltdf$x >= 380 & pltdf$x <= 402 & pltdf$y >= 423 & pltdf$y <= 453),]
d6 = pltdf[which(pltdf$x >= 370 & pltdf$x <= 398 & pltdf$y >= 476 & pltdf$y <= 508),]


##提取区域目标细胞的spot
cltp = c("Sensory Epithelium") ###change
apex = object2_anno[,colnames(object2_anno) %in% rownames(d1)[d1$region_cluster == cltp]]
apex_se_spots <- colnames(apex)

middle = object2_anno[,colnames(object2_anno) %in% rownames(d2)[d2$region_cluster == cltp]]
middle_se_spots <- colnames(middle)

basal = object2_anno[,colnames(object2_anno) %in% rownames(d3)[d3$region_cluster == cltp]]
basal_se_spots <- colnames(basal)


object2_anno$region_annotation <- object2_anno$region_cluster
##对目标spot进行命名
# 转换为字符型
object2_anno$region_annotation <- as.character(object2_anno$region_annotation)
# 赋值
object2_anno$region_annotation[colnames(object2_anno) %in% apex_se_spots] <- "Apex_SE"
object2_anno$region_annotation[colnames(object2_anno) %in% middle_se_spots] <- "Middle_SE"
object2_anno$region_annotation[colnames(object2_anno) %in% basal_se_spots] <- "Basal_SE"

##提取区域目标细胞的spot
cltp = c("Neuron") ###change
apex = object2_anno[,colnames(object2_anno) %in% rownames(d4)[d4$region_cluster == cltp]]
apex_sgn_spots <- colnames(apex)
middle = object2_anno[,colnames(object2_anno) %in% rownames(d5)[d5$region_cluster == cltp]]
middle_sgn_spots <- colnames(middle)
basal = object2_anno[,colnames(object2_anno) %in% rownames(d6)[d6$region_cluster == cltp]]
basal_sgn_spots <- colnames(basal)

object2_anno$region_annotation[colnames(object2_anno) %in% apex_sgn_spots] <- "Apex_Neuron"
object2_anno$region_annotation[colnames(object2_anno) %in% middle_sgn_spots] <- "Middle_Neuron"
object2_anno$region_annotation[colnames(object2_anno) %in% basal_sgn_spots] <- "Basal_Neuron"
#再转回因子型
object2_anno$region_annotation <- as.factor(object2_anno$region_annotation)

Idents(object2_anno) <- object2_anno$region_annotation
levels(object2_anno)
object2_anno <- subset(object2_anno,idents = c("Outer structure",
                                               "Erythrocytes","Macrophages","Endothelial cell","Neutrophils",
                                               "Fibrocytes","SLg_Fibrocyte","SLb_Fibrocyte",                    
                                               "Apex_Neuron","Middle_Neuron","Basal_Neuron",
                                               "Schwann cell",
                                               "Reissner's membrane",             
                                               "Stria",
                                               "Outer sulcus cell",
                                               "Tympanic border cell",
                                               "Apex_SE","Middle_SE","Basal_SE"))
Idents(object2_anno) <- factor(Idents(object2_anno), levels =c("Outer structure",
                                                               "Erythrocytes","Macrophages","Endothelial cell","Neutrophils",
                                                               "Fibrocytes","SLg_Fibrocyte","SLb_Fibrocyte",                     
                                                               "Apex_Neuron","Middle_Neuron","Basal_Neuron",
                                                               "Schwann cell",
                                                               "Reissner's membrane",             
                                                               "Stria",
                                                               "Outer sulcus cell",
                                                               "Tympanic border cell",
                                                               "Apex_SE","Middle_SE","Basal_SE"))

object2_anno$region_annotation <- Idents(object2_anno)

pltdf = data.frame(object2_anno@meta.data, x=object2_anno@images[["sample1"]]@coordinates[["col"]],
                   y=900-object2_anno@images[["sample1"]]@coordinates[["row"]])

object2_color <- c("#F5F5F5",
                   "#F5F5F5","#F5F5F5","#F5F5F5","#F5F5F5",
                   "#F5F5F5","#F5F5F5","#F5F5F5",#"Fibrocytes","SL_Fibrocytes","LW_fibrocyte"
                   "#A23EA5", "#2E86C1","#1ABC9C",#"Neuron",
                   "#F5F5F5",#"Schwann cells"
                   "#F5F5F5",#"Reissner membrane"
                   "#F5F5F5",#stria
                   "#F5F5F5",
                   "#F5F5F5",
                   "#F5F5F5","#F5F5F5","#F5F5F5")
ggplot(pltdf, 
       aes(x, y, color=region_annotation)) + 
  geom_point( size=1) + 
  scale_color_manual(values = object2_color) +
  theme_bw() +
  theme_light(base_size = 15) +
  labs(title = "") +
  theme(plot.title = element_text(hjust = 0.5))

ggsave("E17.5_SCT_L_SGN_region_select_result.pdf",width = 8.3,height = 6.8)

saveRDS(object2_anno,"E17_SCT_L_anno.Rds")
remove(object2_anno)



############################### For P8 Up ###################
Idents(object3_anno) <- object3_anno$celltype
levels(object3_anno)

ids <- c("Outer structure","Neutrophils","Osteoblasts",
         "SLg_Fibrocyte","SLb_Fibrocyte",                    
         "Neuron","Schwann cell",
         "Reissner's membrane",             
         "Stria","Spiral prominence",
         "Outer sulcus cell",
         "Interdental cell",
         "Sensory Epithelium","Sensory Epithelium","Sensory Epithelium",
         "Tympanic border cell",             
         "Sensory Epithelium",                
         "Sensory Epithelium")
names(ids) <- levels(object3_anno)
object3_anno <- RenameIdents(object3_anno,ids)
object3_anno$region_cluster <- Idents(object3_anno)
pltdf = data.frame(object3_anno@meta.data, x=object3_anno@images[["sample1"]]@coordinates[["col"]],
                   y=900-object3_anno@images[["sample1"]]@coordinates[["row"]])

object3_color <- c("#F7ECFD", "#F5E8FD", "#F3E3FD",
                   "#87C986","#33A02C",
                   "#A23EA5",#"Neuron",
                   "#0C2C84",#"Schwann cells"
                   "#E95C59",#"Reissner membrane"
                   "#B49F49",#stia
                   "#00B9DB",#Root cell,Spindle cells
                   "#DB4C6C",
                   "#FFE4B5",
                   "#CD2027",
                   "#FD8D3C"
                   )
ggplot(pltdf, 
       aes(x, y, color=region_cluster)) + 
  geom_point( size=0.6) + 
  scale_color_manual(values = object3_color) +
  theme_bw() +
  theme_light(base_size = 15) +
  labs(title = "") +
  theme(plot.title = element_text(hjust = 0.5))

p = ggplot(pltdf, 
           aes(x, y, color=region_cluster)) + 
  geom_point( size=0.6) + 
  scale_color_manual(values = object3_color) +
  theme_bw() +
  theme_light(base_size = 15) +
  labs(title = "") +
  theme(plot.title = element_text(hjust = 0.5))

##Sensory Epithelium

p+geom_rect(aes(xmin=430,xmax=490,ymin=450,ymax=640),  #Apex
            fill = NA,linewidth = 1,color="red") +
  geom_rect(aes(xmin=500,xmax=565,ymin=450,ymax=525), #Middle
            fill = NA,linewidth = 1,color="blue") +
  geom_rect(aes(xmin=500,xmax=550,ymin=600,ymax=675), #Basal
            fill = NA,linewidth = 1,color="black")
ggsave("P8_SCT_Up_Sensory Epithelium_region_select.pdf",width = 10,height = 8)
##SGN分区
p + geom_rect(aes(xmin=450,xmax=500,ymin=530,ymax=590), #S1
              fill = NA,linewidth = 1,color="red") +
  #geom_rect(aes(xmin=380,xmax=402,ymin=423,ymax=453), #S2
            #fill = NA,linewidth = 1,color="blue") +
  geom_rect(aes(xmin=525,xmax=575,ymin=530,ymax=610), #S3
            fill = NA,linewidth = 1,color="black")
ggsave("P8_SCT_Up_SGN_region_select.pdf",width = 10,height = 8)

#数据提取
#for SE
d1 = pltdf[which(pltdf$x >= 430 & pltdf$x <= 490 & pltdf$y >= 450 & pltdf$y <= 640),]
d2 = pltdf[which(pltdf$x >= 500 & pltdf$x <= 565 & pltdf$y >= 450 & pltdf$y <= 525),]
d3 = pltdf[which(pltdf$x >= 500 & pltdf$x <= 550 & pltdf$y >= 600 & pltdf$y <= 675),]
#for Neuron
d4 = pltdf[which(pltdf$x >= 450 & pltdf$x <= 500 & pltdf$y >= 530 & pltdf$y <= 590),]
#d5 = pltdf[which(pltdf$x >= 380 & pltdf$x <= 402 & pltdf$y >= 423 & pltdf$y <= 453),]
d6 = pltdf[which(pltdf$x >= 525 & pltdf$x <= 575 & pltdf$y >= 530 & pltdf$y <= 610),]


##提取区域目标细胞的spot
cltp = c("Sensory Epithelium") ###change
apex = object3_anno[,colnames(object3_anno) %in% rownames(d1)[d1$region_cluster == cltp]]
apex_se_spots <- colnames(apex)

middle = object3_anno[,colnames(object3_anno) %in% rownames(d2)[d2$region_cluster == cltp]]
middle_se_spots <- colnames(middle)

basal = object3_anno[,colnames(object3_anno) %in% rownames(d3)[d3$region_cluster == cltp]]
basal_se_spots <- colnames(basal)


object3_anno$region_annotation <- object3_anno$region_cluster
##对目标spot进行命名
# 转换为字符型
object3_anno$region_annotation <- as.character(object3_anno$region_annotation)
# 赋值
object3_anno$region_annotation[colnames(object3_anno) %in% apex_se_spots] <- "Apex_SE"
object3_anno$region_annotation[colnames(object3_anno) %in% middle_se_spots] <- "Middle_SE"
object3_anno$region_annotation[colnames(object3_anno) %in% basal_se_spots] <- "Basal_SE"

##提取区域目标细胞的spot
cltp = c("Neuron") ###change
apex = object3_anno[,colnames(object3_anno) %in% rownames(d4)[d4$region_cluster == cltp]]
apex_sgn_spots <- colnames(apex)
#middle = object3_anno[,colnames(object3_anno) %in% rownames(d5)[d5$region_cluster == cltp]]
#middle_sgn_spots <- colnames(middle)
basal = object3_anno[,colnames(object3_anno) %in% rownames(d6)[d6$region_cluster == cltp]]
basal_sgn_spots <- colnames(basal)

object3_anno$region_annotation[colnames(object3_anno) %in% apex_sgn_spots] <- "Apex_Neuron"
#object3_anno$region_annotation[colnames(object3_anno) %in% middle_sgn_spots] <- "Middle_Neuron"
object3_anno$region_annotation[colnames(object3_anno) %in% basal_sgn_spots] <- "Basal_Neuron"
#再转回因子型
object3_anno$region_annotation <- as.factor(object3_anno$region_annotation)

Idents(object3_anno) <- object3_anno$region_annotation
levels(object3_anno)
object3_anno <- subset(object3_anno,idents = c("Outer structure","Neutrophils","Osteoblasts",
                                               "SLg_Fibrocyte","SLb_Fibrocyte",                     
                                               "Apex_Neuron","Basal_Neuron",
                                               "Schwann cell",
                                               "Reissner's membrane",             
                                               "Stria","Spiral prominence",
                                               "Outer sulcus cell","Interdental cell","Tympanic border cell",
                                               "Apex_SE","Middle_SE","Basal_SE"))
Idents(object3_anno) <- factor(Idents(object3_anno), levels =c("Outer structure","Neutrophils","Osteoblasts",
                                                               "SLg_Fibrocyte","SLb_Fibrocyte",                     
                                                               "Apex_Neuron","Basal_Neuron",
                                                               "Schwann cell",
                                                               "Reissner's membrane",             
                                                               "Stria","Spiral prominence",
                                                               "Outer sulcus cell","Interdental cell","Tympanic border cell",
                                                               "Apex_SE","Middle_SE","Basal_SE"))

object3_anno$region_annotation <- Idents(object3_anno)

pltdf = data.frame(object3_anno@meta.data, x=object3_anno@images[["sample1"]]@coordinates[["col"]],
                   y=900-object3_anno@images[["sample1"]]@coordinates[["row"]])

object3_color <- c("#F5F5F5",
                   "#F5F5F5","#F5F5F5",
                   "#F5F5F5","#F5F5F5",#"Fibrocytes","SL_Fibrocytes","LW_fibrocyte"
                   "#A23EA5","#1ABC9C",#"Neuron",
                   "#F5F5F5",#"Schwann cells"
                   "#F5F5F5",#"Reissner membrane"
                   "#F5F5F5",#stria
                   "#F5F5F5",
                   "#F5F5F5",
                   "#F5F5F5",
                   "#F5F5F5",
                   "#F5F5F5","#F5F5F5","#F5F5F5")
ggplot(pltdf, 
       aes(x, y, color=region_annotation)) + 
  geom_point( size=1) + 
  scale_color_manual(values = object3_color) +
  theme_bw() +
  theme_light(base_size = 15) +
  labs(title = "") +
  theme(plot.title = element_text(hjust = 0.5))

ggsave("P8_SCT_Up_SGN_region_select_result.pdf",width = 10,height = 8)

saveRDS(object3_anno,"P8_SCT_Up_anno.Rds")
remove(object3_anno)
################################# P8 down #########################

Idents(object4_anno) <- object4_anno$celltype
levels(object4_anno)
ids <- c("Outer structure",
         "Neutrophils","Chondrocyte","Pericytes","Osteoblasts",
         "SLg_Fibrocyte","SLb_Fibrocyte",                    
         "Neuron","Schwann cell",
         "Reissner's membrane",             
         "Stria","Spiral prominence","Outer sulcus cell","Interdental cell",
         "Sensory Epithelium","Sensory Epithelium",
         "Sensory Epithelium",
         "Tympanic border cell",
         "Sensory Epithelium",                
         "Sensory Epithelium")
names(ids) <- levels(object4_anno)
object4_anno <- RenameIdents(object4_anno,ids)
object4_anno$region_cluster <- Idents(object4_anno)
pltdf = data.frame(object4_anno@meta.data, x=object4_anno@images[["sample1"]]@coordinates[["col"]],
                   y=900-object4_anno@images[["sample1"]]@coordinates[["row"]])

object4_color <- c("#F7ECFD", "#F5E8FD", "#F3E3FD","#F1DFFC","#EFDAFC",
                   "#87C986","#33A02C",
                   "#A23EA5",#"Neuron",
                   "#0C2C84",#"Schwann cells"
                   "#E95C59",#"Reissner membrane"
                   "#B49F49",#stria
                   "#00B9DB",#Root cell,Spindle cells
                   "#DB4C6C","#FFE4B5",
                   "#CD2027","#FD8D3C"
                   )
ggplot(pltdf, 
       aes(x, y, color=region_cluster)) + 
  geom_point( size=0.6) + 
  scale_color_manual(values = object4_color) +
  theme_bw() +
  theme_light(base_size = 15) +
  labs(title = "") +
  theme(plot.title = element_text(hjust = 0.5))

p = ggplot(pltdf, 
           aes(x, y, color=region_cluster)) + 
  geom_point( size=0.6) + 
  scale_color_manual(values = object4_color) +
  theme_bw() +
  theme_light(base_size = 15) +
  labs(title = "") +
  theme(plot.title = element_text(hjust = 0.5))

##Sensory Epithelium

p+geom_rect(aes(xmin=430,xmax=525,ymin=90,ymax=260),  #Apex
            fill = NA,linewidth = 1,color="red") +
  geom_rect(aes(xmin=525,xmax=600,ymin=100,ymax=150), #Middle
            fill = NA,linewidth = 1,color="blue") +
  geom_rect(aes(xmin=525,xmax=600,ymin=250,ymax=300), #Basal
            fill = NA,linewidth = 1,color="black")
ggsave("P8_SCT_Down_Sensory Epithelium_region_select.pdf",width = 10,height = 8)
##SGN分区
p + geom_rect(aes(xmin=500,xmax=550,ymin=170,ymax=220), #S1
              fill = NA,linewidth = 1,color="red") +
  geom_rect(aes(xmin=560,xmax=610,ymin=140,ymax=180), #S2
  fill = NA,linewidth = 1,color="blue") +
  geom_rect(aes(xmin=580,xmax=620,ymin=210,ymax=265), #S3
            fill = NA,linewidth = 1,color="black")
ggsave("P8_SCT_Down_SGN_region_select.pdf",width = 10.3,height = 8.3)

#数据提取
#for SE
d1 = pltdf[which(pltdf$x >= 430 & pltdf$x <= 525 & pltdf$y >= 90 & pltdf$y <= 260),]
d2 = pltdf[which(pltdf$x >= 525 & pltdf$x <= 600 & pltdf$y >= 100 & pltdf$y <= 150),]
d3 = pltdf[which(pltdf$x >= 525 & pltdf$x <= 600 & pltdf$y >= 250 & pltdf$y <= 300),]
#for Neuron
d4 = pltdf[which(pltdf$x >= 500 & pltdf$x <= 550 & pltdf$y >= 170 & pltdf$y <= 220),]
d5 = pltdf[which(pltdf$x >= 560 & pltdf$x <= 610 & pltdf$y >= 140 & pltdf$y <= 180),]
d6 = pltdf[which(pltdf$x >= 580 & pltdf$x <= 620 & pltdf$y >= 210 & pltdf$y <= 265),]


##提取区域目标细胞的spot
cltp = c("Sensory Epithelium") ###change
apex = object4_anno[,colnames(object4_anno) %in% rownames(d1)[d1$region_cluster == cltp]]
apex_se_spots <- colnames(apex)

middle = object4_anno[,colnames(object4_anno) %in% rownames(d2)[d2$region_cluster == cltp]]
middle_se_spots <- colnames(middle)

basal = object4_anno[,colnames(object4_anno) %in% rownames(d3)[d3$region_cluster == cltp]]
basal_se_spots <- colnames(basal)


object4_anno$region_annotation <- object4_anno$region_cluster
##对目标spot进行命名
# 转换为字符型
object4_anno$region_annotation <- as.character(object4_anno$region_annotation)
# 赋值
object4_anno$region_annotation[colnames(object4_anno) %in% apex_se_spots] <- "Apex_SE"
object4_anno$region_annotation[colnames(object4_anno) %in% middle_se_spots] <- "Middle_SE"
object4_anno$region_annotation[colnames(object4_anno) %in% basal_se_spots] <- "Basal_SE"

##提取区域目标细胞的spot
cltp = c("Neuron") ###change
apex = object4_anno[,colnames(object4_anno) %in% rownames(d4)[d4$region_cluster == cltp]]
apex_sgn_spots <- colnames(apex)
middle = object4_anno[,colnames(object4_anno) %in% rownames(d5)[d5$region_cluster == cltp]]
middle_sgn_spots <- colnames(middle)
basal = object4_anno[,colnames(object4_anno) %in% rownames(d6)[d6$region_cluster == cltp]]
basal_sgn_spots <- colnames(basal)

object4_anno$region_annotation[colnames(object4_anno) %in% apex_sgn_spots] <- "Apex_Neuron"
object4_anno$region_annotation[colnames(object4_anno) %in% middle_sgn_spots] <- "Middle_Neuron"
object4_anno$region_annotation[colnames(object4_anno) %in% basal_sgn_spots] <- "Basal_Neuron"
#再转回因子型
object4_anno$region_annotation <- as.factor(object4_anno$region_annotation)

Idents(object4_anno) <- object4_anno$region_annotation
levels(object4_anno)
object4_anno <- subset(object4_anno,idents = c("Outer structure",
                                               "Neutrophils","Chondrocyte","Pericytes","Osteoblasts",
                                               "SLg_Fibrocyte","SLb_Fibrocyte",                       
                                               "Apex_Neuron","Middle_Neuron","Basal_Neuron",
                                               "Schwann cell",
                                               "Reissner's membrane",             
                                               "Stria","Spiral prominence","Outer sulcus cell","Interdental cell",
                                               "Tympanic border cell",
                                               "Apex_SE","Middle_SE","Basal_SE"))
Idents(object4_anno) <- factor(Idents(object4_anno), levels =c("Outer structure",
                                                               "Neutrophils","Chondrocyte","Pericytes","Osteoblasts",
                                                               "SLg_Fibrocyte","SLb_Fibrocyte",                       
                                                               "Apex_Neuron","Middle_Neuron","Basal_Neuron",
                                                               "Schwann cell",
                                                               "Reissner's membrane",             
                                                               "Stria","Spiral prominence","Outer sulcus cell","Interdental cell",
                                                               "Tympanic border cell",
                                                               "Apex_SE","Middle_SE","Basal_SE"))

object4_anno$region_annotation <- Idents(object4_anno)

pltdf = data.frame(object4_anno@meta.data, x=object4_anno@images[["sample1"]]@coordinates[["col"]],
                   y=900-object4_anno@images[["sample1"]]@coordinates[["row"]])

object4_color <- c("#F5F5F5",
                   "#F5F5F5","#F5F5F5","#F5F5F5","#F5F5F5",
                   "#F5F5F5","#F5F5F5",#"SL_Fibrocytes","LW_fibrocyte"
                   "#A23EA5", "#2E86C1","#1ABC9C",#"Neuron",
                   "#F5F5F5",#"Schwann cells"
                   "#F5F5F5",#"Reissner membrane"
                   "#F5F5F5",#stria
                   "#F5F5F5",#SP
                   "#F5F5F5",#OSC
                   "#F5F5F5",#IdC
                   "#F5F5F5",#TBC
                   "#F5F5F5","#F5F5F5","#F5F5F5")
ggplot(pltdf, 
       aes(x, y, color=region_annotation)) + 
  geom_point( size=1) + 
  scale_color_manual(values = object4_color) +
  theme_bw() +
  theme_light(base_size = 15) +
  labs(title = "") +
  theme(plot.title = element_text(hjust = 0.5))

ggsave("P8_SCT_Down_SGN_region_select_result.pdf",width = 11.5,height = 9.5)

saveRDS(object4_anno,"P8_SCT_Down_anno.Rds")
remove(object4_anno)
####################################### Adult right ######################
Idents(object6_anno) <- object6_anno$celltype

levels(object6_anno)

ids <- c("Outer structure",
         "Pericytes","Neutrophils","Osteoblasts","Endothelial cell", 
         "SLg_Fibrocyte","SLb_Fibrocyte",                    
         "Neuron","Schwann cell",
         "Reissner's membrane",             
         "Stria","Spiral prominence",
         "Interdental cell","Sensory Epithelium",
         "Tympanic border cell",
         "Sensory Epithelium",                
         "Sensory Epithelium")
names(ids) <- levels(object6_anno)
object6_anno <- RenameIdents(object6_anno,ids)
object6_anno$region_cluster <- Idents(object6_anno)
pltdf = data.frame(object6_anno@meta.data, x=object6_anno@images[["sample1"]]@coordinates[["col"]],
                   y=900-object6_anno@images[["sample1"]]@coordinates[["row"]])

object6_color <- c("#F7ECFD", "#F5E8FD", "#F3E3FD",
                   "#F1DFFC","#EFDAFC",
                   "#C5DEBA","#33A02C",
                   "#A23EA5",#"Neuron",
                   "#0C2C84",#"Schwann cells"
                   "#E95C59",#"Reissner membrane"
                   "#B49F49",#stia
                   "#00B9DB",#Root cell,Spindle cells
                   "#DB4C6C",
                   "#CD2027",#OHC,IHC
                   "#FFE4B5"
)
ggplot(pltdf, 
       aes(x, y, color=region_cluster)) + 
  geom_point( size=0.6) + 
  scale_color_manual(values = object6_color) +
  theme_bw() +
  theme_light(base_size = 15) +
  labs(title = "") +
  theme(plot.title = element_text(hjust = 0.5))

p = ggplot(pltdf, 
           aes(x, y, color=region_cluster)) + 
  geom_point( size=0.6) + 
  scale_color_manual(values = object6_color) +
  theme_bw() +
  theme_light(base_size = 15) +
  labs(title = "") +
  theme(plot.title = element_text(hjust = 0.5))

##Sensory Epithelium

p+geom_rect(aes(xmin=590,xmax=760,ymin=350,ymax=400),  #Apex
            fill = NA,linewidth = 1,color="red") +
  geom_rect(aes(xmin=575,xmax=650,ymin=260,ymax=320), #Middle
            fill = NA,linewidth = 1,color="blue") +
  geom_rect(aes(xmin=750,xmax=810,ymin=290,ymax=325), #Basal
            fill = NA,linewidth = 1,color="black")
ggsave("Adult_SCT_Right_Sensory Epithelium_region_select.pdf",width = 10,height = 8)
##SGN分区
p + geom_rect(aes(xmin=650,xmax=720,ymin=325,ymax=370), #S1
              fill = NA,linewidth = 1,color="red") +
  geom_rect(aes(xmin=625,xmax=690,ymin=260,ymax=310), #S2
            fill = NA,linewidth = 1,color="blue") +
  geom_rect(aes(xmin=700,xmax=760,ymin=260,ymax=305), #S3
            fill = NA,linewidth = 1,color="black")
ggsave("Adult_SCT_Right_SGN_region_select.pdf",width = 10,height = 7.1)

#数据提取
#for SE
d1 = pltdf[which(pltdf$x >= 590 & pltdf$x <= 760 & pltdf$y >= 350 & pltdf$y <= 400),]
d2 = pltdf[which(pltdf$x >= 575 & pltdf$x <= 650 & pltdf$y >= 260 & pltdf$y <= 320),]
d3 = pltdf[which(pltdf$x >= 750 & pltdf$x <= 810 & pltdf$y >= 290 & pltdf$y <= 325),]
#for Neuron
d4 = pltdf[which(pltdf$x >= 650 & pltdf$x <= 720 & pltdf$y >= 325 & pltdf$y <= 370),]
d5 = pltdf[which(pltdf$x >= 625 & pltdf$x <= 690 & pltdf$y >= 260 & pltdf$y <= 310),]
d6 = pltdf[which(pltdf$x >= 700 & pltdf$x <= 760 & pltdf$y >= 260 & pltdf$y <= 305),]


##提取区域目标细胞的spot
cltp = c("Sensory Epithelium") ###change
apex = object6_anno[,colnames(object6_anno) %in% rownames(d1)[d1$region_cluster == cltp]]
apex_se_spots <- colnames(apex)

middle = object6_anno[,colnames(object6_anno) %in% rownames(d2)[d2$region_cluster == cltp]]
middle_se_spots <- colnames(middle)

basal = object6_anno[,colnames(object6_anno) %in% rownames(d3)[d3$region_cluster == cltp]]
basal_se_spots <- colnames(basal)


object6_anno$region_annotation <- object6_anno$region_cluster
##对目标spot进行命名
# 转换为字符型
object6_anno$region_annotation <- as.character(object6_anno$region_annotation)
# 赋值
object6_anno$region_annotation[colnames(object6_anno) %in% apex_se_spots] <- "Apex_SE"
object6_anno$region_annotation[colnames(object6_anno) %in% middle_se_spots] <- "Middle_SE"
object6_anno$region_annotation[colnames(object6_anno) %in% basal_se_spots] <- "Basal_SE"

##提取区域目标细胞的spot
cltp = c("Neuron") ###change
apex = object6_anno[,colnames(object6_anno) %in% rownames(d4)[d4$region_cluster == cltp]]
apex_sgn_spots <- colnames(apex)
middle = object6_anno[,colnames(object6_anno) %in% rownames(d5)[d5$region_cluster == cltp]]
middle_sgn_spots <- colnames(middle)
basal = object6_anno[,colnames(object6_anno) %in% rownames(d6)[d6$region_cluster == cltp]]
basal_sgn_spots <- colnames(basal)

object6_anno$region_annotation[colnames(object6_anno) %in% apex_sgn_spots] <- "Apex_Neuron"
object6_anno$region_annotation[colnames(object6_anno) %in% middle_sgn_spots] <- "Middle_Neuron"
object6_anno$region_annotation[colnames(object6_anno) %in% basal_sgn_spots] <- "Basal_Neuron"
#再转回因子型
object6_anno$region_annotation <- as.factor(object6_anno$region_annotation)

Idents(object6_anno) <- object6_anno$region_annotation
levels(object6_anno)
object6_anno <- subset(object6_anno,idents = c("Outer structure",
                                               "Pericytes","Neutrophils","Osteoblasts","Endothelial cell", 
                                               "SLg_Fibrocyte","SLb_Fibrocyte",                       
                                               "Apex_Neuron","Middle_Neuron","Basal_Neuron",
                                               "Schwann cell",
                                               "Reissner's membrane",             
                                               "Stria","Spiral prominence",
                                               "Interdental cell",
                                               "Tympanic border cell",
                                               "Apex_SE","Middle_SE","Basal_SE"))
Idents(object6_anno) <- factor(Idents(object6_anno), levels =c("Outer structure",
                                                               "Pericytes","Neutrophils","Osteoblasts","Endothelial cell", 
                                                               "SLg_Fibrocyte","SLb_Fibrocyte",                       
                                                               "Apex_Neuron","Middle_Neuron","Basal_Neuron",
                                                               "Schwann cell",
                                                               "Reissner's membrane",             
                                                               "Stria","Spiral prominence",
                                                               "Interdental cell",
                                                               "Tympanic border cell",
                                                               "Apex_SE","Middle_SE","Basal_SE"))

object6_anno$region_annotation <- Idents(object6_anno)

pltdf = data.frame(object6_anno@meta.data, x=object6_anno@images[["sample1"]]@coordinates[["col"]],
                   y=900-object6_anno@images[["sample1"]]@coordinates[["row"]])

object6_color <- c("#F5F5F5",
                   "#F5F5F5","#F5F5F5","#F5F5F5","#F5F5F5",
                   "#F5F5F5","#F5F5F5",#"SL_Fibrocytes","LW_fibrocyte"
                   "#A23EA5", "#2E86C1","#1ABC9C",#"Neuron",
                   "#F5F5F5",#"Schwann cells"
                   "#F5F5F5",#"Reissner membrane"
                   "#F5F5F5",#stria
                   "#F5F5F5",#SP
                   "#F5F5F5",#IdC
                   "#F5F5F5",#
                   "#F5F5F5","#F5F5F5","#F5F5F5")
ggplot(pltdf, 
       aes(x, y, color=region_annotation)) + 
  geom_point( size=1) + 
  scale_color_manual(values = object6_color) +
  theme_bw() +
  theme_light(base_size = 15) +
  labs(title = "") +
  theme(plot.title = element_text(hjust = 0.5))

ggsave("Adult_SCT_Right_SGN_region_select_result.pdf",width = 11.7,height = 8.3)

saveRDS(object6_anno,"Adult_SCT_Right_anno.Rds")
remove(object6_anno)

####################不同分区SGN数量统计###################

object1_region <- data.frame(region=object1_anno$region_annotation,
                             maturation_score=object1_anno$Pseudotime)
object1_region$cell_ID <- rownames(object1_region)
object1_region$sample <- "E17.5_Right"
object1_region$group <- "E17.5"
object1_region_SGN <- object1_region[object1_region$region %in% c("Apex_Neuron","Middle_Neuron","Basal_Neuron"),]

object2_region <- data.frame(region=object2_anno$region_annotation,
                             maturation_score=object2_anno$Pseudotime)
object2_region$cell_ID <- rownames(object2_region)
object2_region$sample <- "E17.5_Left"
object2_region$group <- "E17.5"
object2_region_SGN <- object2_region[object2_region$region %in% c("Apex_Neuron","Middle_Neuron","Basal_Neuron"),]

object3_region <- data.frame(region=object3_anno$region_annotation,
                             maturation_score=object3_anno$Pseudotime)
object3_region$cell_ID <- rownames(object3_region)
object3_region$sample <- "P8_Up"
object3_region$group <- "P8"
object3_region_SGN <- object3_region[object3_region$region %in% c("Apex_Neuron","Middle_Neuron","Basal_Neuron"),]

object4_region <- data.frame(region=object4_anno$region_annotation,
                             maturation_score=object4_anno$Pseudotime)
object4_region$cell_ID <- rownames(object4_region)
object4_region$sample <- "P8_Down"
object4_region$group <- "P8"
object4_region_SGN <- object4_region[object4_region$region %in% c("Apex_Neuron","Middle_Neuron","Basal_Neuron"),]

object6_region <- data.frame(region=object6_anno$region_annotation,
                             maturation_score=object6_anno$Pseudotime)
object6_region$cell_ID <- rownames(object6_region)
object6_region$sample <- "Adult_Right"
object6_region$group <- "Adult"
object6_region_SGN <- object6_region[object6_region$region %in% c("Apex_Neuron","Middle_Neuron","Basal_Neuron"),]

SGN_region_all <- rbind(object1_region_SGN,object2_region_SGN,object3_region_SGN,object4_region_SGN,object6_region_SGN)
SGN_region_all_new <- data.frame(sample=SGN_region_all$sample,
                                 group=SGN_region_all$group,
                             cell=SGN_region_all$cell_ID,
                             region=SGN_region_all$region,
                             maturation_score=SGN_region_all$maturation_score)
write.csv(SGN_region_all_new,"SCT_all_Neuron_region_cell_count_mature_score.csv")
##统计不同分区SGN数量
plot_data <- SGN_region_all_new %>%
  group_by(sample, region) %>%
  summarise(count = n(), .groups = "drop")
plot_data <- plot_data %>%
  mutate(sample = factor(sample, levels = c("E17.5_Right", "E17.5_Left",
                                            "P8_Up", "P8_Down",
                                            "Adult_Right")))
region_colors_highlight <- c("Apex_Neuron"="#A23EA5",
                             "Middle_Neuron"="#2E86C1",
                             "Basal_Neuron"="#1ABC9C")
ggplot(plot_data, aes(x = sample, y = count, fill = region)) +
  geom_col(width=0.8,position = position_dodge(width = 0.8,preserve = "single")) +  # 分组柱状
  scale_fill_manual(values = region_colors_highlight) +
  labs(
    y = "Cell Number",
    fill = "Region"
  ) +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        axis.title.x = element_blank(),
        axis.line = element_line(linewidth = 1),
        axis.ticks = element_line(linewidth = 1),
        axis.text = element_text(size = 10),
    axis.text.x = element_text(angle = 45, hjust = 1)  # 样本名倾斜防重叠
  )
ggsave("SCT_allsample_Neuron_region_cell_count.pdf",width = 6,height = 4)
##统计不同时期不同分区SGN数量
plot_data <- SGN_region_all_new %>%
  group_by(group, region) %>%
  summarise(count = n(), .groups = "drop")
plot_data <- plot_data %>%
  mutate(group = factor(group, levels = c("E17.5", 
                                            "P8",
                                            "Adult")))
region_colors_highlight <- c("Apex_Neuron"="#A23EA5",
                             "Middle_Neuron"="#2E86C1",
                             "Basal_Neuron"="#1ABC9C")
ggplot(plot_data, aes(x = group, y = count, fill = region)) +
  geom_col(width=0.8,position = position_dodge(width = 0.8,preserve = "single")) +  # 分组柱状
  scale_fill_manual(values = region_colors_highlight) +
  labs(
    y = "Cell Number",
    fill = "Region"
  ) +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        axis.title.x = element_blank(),
        axis.line = element_line(linewidth = 1),
        axis.ticks = element_line(linewidth = 1),
        axis.text = element_text(size = 10),
        axis.text.x = element_text(angle = 45, hjust = 1)  # 样本名倾斜防重叠
  )
ggsave("SCT_allgroup_Neuron_region_cell_count.pdf",width = 6,height = 4)

##统计不同样本不同分区SGN成熟度分数
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(rstatix)
library(patchwork)
library(multcompView)
library(RColorBrewer)
maturation_stats <- SGN_region_all_new %>%
  group_by(sample, region) %>%
  summarise(
    mean_score = mean(maturation_score, na.rm = TRUE),
    sd_score = sd(maturation_score, na.rm = TRUE),
    count = n(),
    .groups = "drop"
  )
maturation_stats <- maturation_stats %>%
  mutate(sample = factor(sample, levels = c("E17.5_Right", "E17.5_Left",
                                            "P8_Up", "P8_Down",
                                            "Adult_Right")))


anova_by_group <- SGN_region_all_new %>%
  group_by(group) %>%
  anova_test(maturation_score ~ region) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance()


region_colors_highlight <- c("Apex_Neuron"="#A23EA5",
                             "Middle_Neuron"="#2E86C1",
                             "Basal_Neuron"="#1ABC9C")
ggplot(maturation_stats, aes(x = sample, y = mean_score, fill = region)) +
  geom_col(width=0.8,position = position_dodge(width = 0.8,preserve = "single")) +  # 分组柱状
  geom_bar(stat = "identity", position = position_dodge(0.9)) +
  geom_errorbar(aes(ymin = mean_score - sd_score, ymax = mean_score + sd_score),
                width = 0.2, position = position_dodge(0.9)) +
  scale_fill_manual(values = region_colors_highlight) +
  labs(
    y = "Maturation score",
    fill = "Region"
  ) +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        axis.title.x = element_blank(),
        axis.line = element_line(linewidth = 1),
        axis.ticks = element_line(linewidth = 1),
        axis.text = element_text(size = 10),
        axis.text.x = element_text(angle = 45, hjust = 1)  # 样本名倾斜防重叠
  )
ggsave("SCT_allsample_Neuron_region_cell_Maturation score.pdf",width = 6,height = 4)