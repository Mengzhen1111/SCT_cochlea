library(Seurat)
library(SeuratObject)
library(Seurat)
library(tidyverse)
library(Matrix)
library(spacexr)
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
#library(clustree)
library(patchwork)
library(RColorBrewer)
library(colorRamp2)
library(ComplexHeatmap)
library(semla)
library(ggforce)
library(ggplot2)
library(cowplot)
library(Nebulosa)

###################################### Data QC #################################

# ===================================  E17.5========================
object1 <- readRDS("E17_SCT_R.Rds")
object1 <- PercentageFeatureSet(object1, "^mt-", col.name = "percent_mito")
saveRDS(object1,"E17_SCT_R.Rds")
semla.data1 <- UpdateSeuratForSemla(object1)
semla.data1 <- LoadImages(semla.data1, verbose = FALSE)
MapFeaturesSummary(semla.data1, features = "percent_mito", ncol = 1,
          subplot_type = "violin",
          pt_alpha = 5,
          pt_size = 1,
          pt_stroke = 0.001,
          crop_area = c(0.55, 0.55, 0.75, 0.71)
)&theme(legend.position = "right") 

genes_to_plot <- c("nCount_Spatial", "nFeature_Spatial", "percent_mito")

plot_list <- lapply(genes_to_plot, function(g){
  MapFeaturesSummary(semla.data1, features = g, ncol = 1,
                     subplot_type = "violin",
                     pt_alpha = 5,
                     pt_size = 1,
                     pt_stroke = 0.001,
                     crop_area = c(0.55, 0.55, 0.75, 0.71)
  )&theme(legend.position = "right") 
})
pdf("SCT_E17.5_R_QC.pdf",width =18,height = 5)
wrap_plots(plot_list, nrow = 1)
dev.off()

object1 <- readRDS("E17_SCT_L.Rds")
object1 <- PercentageFeatureSet(object1, "^mt-", col.name = "percent_mito")
saveRDS(object1,"E17_SCT_L.Rds")

semla.data1 <- UpdateSeuratForSemla(object1)
semla.data1 <- LoadImages(semla.data1, verbose = FALSE)
MapFeaturesSummary(semla.data1, features = "nCount_Spatial", ncol = 1,
                   subplot_type = "violin",
                   pt_alpha = 5,
                   pt_size = 1,
                   pt_stroke = 0.001,
                   crop_area = c(0.33, 0.34, 0.49, 0.53)
)&theme(legend.position = "right") 

genes_to_plot <- c("nCount_Spatial", "nFeature_Spatial", "percent_mito")

plot_list <- lapply(genes_to_plot, function(g){
  MapFeaturesSummary(semla.data1, features = g, ncol = 1,
                     subplot_type = "violin",
                     pt_alpha = 5,
                     pt_size = 1,
                     pt_stroke = 0.001,
                     crop_area = c(0.33, 0.34, 0.49, 0.53)
  )&theme(legend.position = "right") 
})
wrap_plots(plot_list, nrow = 1)
pdf("SCT_E17.5_L_QC.pdf",width =17,height = 5)
wrap_plots(plot_list, nrow = 1)
dev.off()

# ===================================  P8 ========================
object1 <- readRDS("P8_SCT_up_anno.Rds")
object1 <- PercentageFeatureSet(object1, "^mt-", col.name = "percent_mito")
saveRDS(object1,"P8_SCT_up_anno.Rds")

semla.data1 <- UpdateSeuratForSemla(object1)
semla.data1 <- LoadImages(semla.data1, verbose = FALSE)
MapFeaturesSummary(semla.data1, features = "percent_mito", ncol = 1,
                   subplot_type = "violin",
                   pt_alpha = 5,
                   pt_size = 1,
                   pt_stroke = 0.001,
                   crop_area = c(0.4, 0.18, 0.65, 0.45)
)&theme(legend.position = "right")

genes_to_plot <- c("nCount_Spatial", "nFeature_Spatial", "percent_mito")

plot_list <- lapply(genes_to_plot, function(g){
  MapFeaturesSummary(semla.data1, features = g, ncol = 1,
                     subplot_type = "violin",
                     pt_alpha = 5,
                     pt_size = 0.9,
                     pt_stroke = 0.001,
                     crop_area = c(0.4, 0.18, 0.65, 0.45)
  )&theme(legend.position = "right") 
})
pdf("SCT_P8_Up_QC.pdf",width =20,height = 5)
wrap_plots(plot_list, nrow = 1)
dev.off()

object1 <- readRDS("P8_SCT_Down_anno.Rds")
object1 <- PercentageFeatureSet(object1, "^mt-", col.name = "percent_mito")
saveRDS(object1,"P8_SCT_Down_anno.Rds")

semla.data1 <- UpdateSeuratForSemla(object1)
semla.data1 <- LoadImages(semla.data1, verbose = FALSE)
MapFeaturesSummary(semla.data1, features = "nCount_Spatial", ncol = 1,
                   subplot_type = "violin",
                   pt_alpha = 5,
                   pt_size = 0.9,
                   pt_stroke = 0.001,
                   crop_area = c(0.45, 0.55, 0.715, 0.835)
)&theme(legend.position = "right") 

genes_to_plot <- c("nCount_Spatial", "nFeature_Spatial", "percent_mito")

plot_list <- lapply(genes_to_plot, function(g){
  MapFeaturesSummary(semla.data1, features = g, ncol = 1,
                     subplot_type = "violin",
                     pt_alpha = 5,
                     pt_size = 0.9,
                     pt_stroke = 0.001,
                     crop_area = c(0.45, 0.55, 0.715, 0.835)
  )&theme(legend.position = "right") 
})
wrap_plots(plot_list, nrow = 1)
pdf("SCT_P8_Down_QC.pdf",width =20,height = 5)
wrap_plots(plot_list, nrow = 1)
dev.off()

# =================================== Adult ========================
object1 <- readRDS("Adult_SCT_right_anno.Rds")
object1 <- PercentageFeatureSet(object1, "^mt-", col.name = "percent_mito")
saveRDS(object1,"Adult_SCT_right_anno.Rds")

semla.data1 <- UpdateSeuratForSemla(object1)
semla.data1 <- LoadImages(semla.data1, verbose = FALSE)
MapFeaturesSummary(semla.data1, features = "percent_mito", ncol = 1,
                   subplot_type = "violin",
                   pt_alpha = 5,
                   pt_size = 0.9,
                   pt_stroke = 0.001,
                   crop_area = c(0.55, 0.49, 0.85, 0.685)
)&theme(legend.position = "right")

genes_to_plot <- c("nCount_Spatial", "nFeature_Spatial", "percent_mito")

plot_list <- lapply(genes_to_plot, function(g){
  MapFeaturesSummary(semla.data1, features = g, ncol = 1,
                     subplot_type = "violin",
                     pt_alpha = 5,
                     pt_size = 0.8,
                     pt_stroke = 0.001,
                     crop_area = c(0.55, 0.49, 0.85, 0.685)
  )&theme(legend.position = "right") 
})
pdf("SCT_Adult_Right_QC.pdf",width =20,height = 5)
wrap_plots(plot_list, nrow = 1)
dev.off()
