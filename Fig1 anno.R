

####################################  Fig1 related code  #########################

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
library(clustree)
library(patchwork)
library(RColorBrewer)
library(colorRamp2)
library(ComplexHeatmap)
library(semla)
library(ggforce)
library(ggplot2)
library(cowplot)
library(Nebulosa)
source("CreateBmkObject.R")
Rcpp::sourceCpp(code='
#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::export]]
IntegerMatrix asMatrix(NumericVector rp,
                       NumericVector cp,
                       NumericVector z,
                       int nrows,
                       int ncols){
  int k = z.size() ;
  IntegerMatrix  mat(nrows, ncols);
  for (int i = 0; i < k; i++){
      mat(rp[i],cp[i]) = z[i];
  }
  return mat;
}
' )
as_matrix <- function(mat){
  
  row_pos <- mat@i
  col_pos <- findInterval(seq(mat@x)-1,mat@p[-1])
  
  tmp <- asMatrix(rp = row_pos, cp = col_pos, z = mat@x,
                  nrows =  mat@Dim[1], ncols = mat@Dim[2])
  
  row.names(tmp) <- mat@Dimnames[[1]]
  colnames(tmp) <- mat@Dimnames[[2]]
  return(tmp)
}
###################################### Data preprocessing #################################
object1 <- CreateS1000Object(
  matrix_path="L2_E17.5_R",
  png_path="he_roi_small.png",
  spot_radius =1,
  min.cells =5,min.features = 50)
object1 <- NormalizeData(object1, normalization.method = "LogNormalize", scale.factor = 10000)
object1 <- FindVariableFeatures(object1, selection.method = "vst", nfeatures = 2000)
object1 <- ScaleData(object1)
object1 <- RunPCA(object=object1,pc.genes = VariableFeatures(object2))
object1 <- FindNeighbors(object=object1, reduction = "pca",dims = 1:30,verbose = F)
object1 <- FindClusters(object=object1, verbose = F,resolution = 1)
object1 <- RunUMAP(object=object1, reduction = "pca",dims = 1:14)
DimPlot(object1, reduction = "umap", label = T,pt.size=1)+NoLegend()
object2$sample <- "E17.5_R"

object2 <- CreateS1000Object(
  matrix_path="L2_E17.5_L",
  png_path="he_roi_small.png",
  spot_radius =1,
  min.cells =5,min.features = 50)
object2 <- NormalizeData(object2, normalization.method = "LogNormalize", scale.factor = 10000)
object2 <- FindVariableFeatures(object2, selection.method = "vst", nfeatures = 2000)
object2 <- ScaleData(object2)
object2 <- RunPCA(object=object2,pc.genes = VariableFeatures(object2))
object2 <- FindNeighbors(object=object2, reduction = "pca",dims = 1:30,verbose = F)
object2 <- FindClusters(object=object2, verbose = F,resolution = 1)
object2 <- RunUMAP(object=object2, reduction = "pca",dims = 1:14)
DimPlot(object2, reduction = "umap", label = T,pt.size=1)+NoLegend()
object2$sample <- "E17.5_L"
saveRDS(object2,"E17_SCT_L.Rds")

object3 <- CreateS1000Object(
  matrix_path="L2_P8_select_up",
  png_path="he_roi_small.png",
  spot_radius =1,
  min.cells =5,min.features = 50)
object3 <- NormalizeData(object3, normalization.method = "LogNormalize", scale.factor = 10000)
object3 <- FindVariableFeatures(object3, selection.method = "vst", nfeatures = 2000)
object3 <- ScaleData(object3)
object3 <- RunPCA(object=object3,pc.genes = VariableFeatures(object3))
object3 <- FindNeighbors(object=object3, reduction = "pca",dims = 1:30,verbose = F)
object3 <- FindClusters(object=object3, verbose = F,resolution = 1)
object3 <- RunUMAP(object=object3, reduction = "pca",dims = 1:14)
object3$sample <- "P8_Up"
saveRDS(object3,"P8_SCT_Up.Rds")

object4 <- CreateS1000Object(
  matrix_path="L2_P8_select_down",
  png_path="he_roi_small.png",
  spot_radius =1,
  min.cells =5,min.features = 50)
object4 <- NormalizeData(object4, normalization.method = "LogNormalize", scale.factor = 10000)
object4 <- FindVariableFeatures(object4, selection.method = "vst", nfeatures = 2000)
object4 <- ScaleData(object4)
object4 <- RunPCA(object=object4,pc.genes = VariableFeatures(object4))
object4 <- FindNeighbors(object=object4, reduction = "pca",dims = 1:30,verbose = F)
object4 <- FindClusters(object=object4, verbose = F,resolution = 1)
object4 <- RunUMAP(object=object4, reduction = "pca",dims = 1:14)
object4$sample <- "P8_Down"
saveRDS(object4,"P8_SCT_Down.Rds")

object5 <- CreateS1000Object(
  matrix_path="L2_Adult_left",
  png_path="he_roi_small.png",
  spot_radius =1,
  min.cells =1,min.features = 1)
object5 <- NormalizeData(object5, normalization.method = "LogNormalize", scale.factor = 10000)
object5 <- FindVariableFeatures(object5, selection.method = "vst", nfeatures = 2000)
object5 <- ScaleData(object5)
object5 <- RunPCA(object=object5,pc.genes = VariableFeatures(object6))
object5 <- FindNeighbors(object=object5, reduction = "pca",dims = 1:30,verbose = F)
object5 <- FindClusters(object=object5, verbose = F,resolution = 1)
object5 <- RunUMAP(object=object5, reduction = "pca",dims = 1:14)
object5$sample <- "Adult_right"
saveRDS(object5,"Adult_SCT_Down.Rds")

object6 <- CreateS1000Object(
  matrix_path="L2_Adult_right",
  png_path="he_roi_small.png",
  spot_radius =1,
  min.cells =5,min.features = 50)
object6 <- NormalizeData(object6, normalization.method = "LogNormalize", scale.factor = 10000)
object6 <- FindVariableFeatures(object6, selection.method = "vst", nfeatures = 2000)
object6 <- ScaleData(object6)
object6 <- RunPCA(object=object6,pc.genes = VariableFeatures(object6))
object6 <- FindNeighbors(object=object6, reduction = "pca",dims = 1:30,verbose = F)
object6 <- FindClusters(object=object6, verbose = F,resolution = 1)
object6 <- RunUMAP(object=object6, reduction = "pca",dims = 1:14)
object6$sample <- "Adult_right"
saveRDS(object6,"Adult_SCT_right.Rds")

###################################### Data annotation #################################

#参考marker list
pnas_marker <- read_xlsx("PNAS_top50_gene.xlsx")
HS_marker <- read_xlsx("P1_markers.xlsx")
E16_marker <- read_xlsx("E16_markers.xlsx")
sc_top100 <- read.csv("mmtop100_cluster_ident.csv")

###################################### E17.5 annotation #################################

###################################### E17.5_R annotation #################################

celltype_1 <- colnames(pnas_marker)
for (i in celltype_1){
  fe <- pnas_marker[i]
  p1 <- DotPlot(object1, features = fe)+ RotatedAxis()+theme(legend.position = "none")+labs(title = i)
  ggsave(paste( i, "_P8.png", sep=""),width = 12,height = 10)
}
for (i in colnames(E16_marker)){
  fe <- E16_marker[i]
  p1 <- DotPlot(object1, features = fe)+ RotatedAxis()+theme(legend.position = "none")+labs(title = i)
  ggsave(paste( i, "_E16.png", sep=""),width = 12,height = 10)
}


for (i in colnames(sc_top100)){
  fe <- sc_top100[i]
  p1 <- DotPlot(object1, features = fe)+ RotatedAxis()+theme(legend.position = "none")+labs(title = i)
  ggsave(paste( i, "_adult.png", sep=""),width = 12,height = 10)
}
for (i in colnames(HS_marker)){
  fe <- HS_marker[i]
  p1 <- DotPlot(object1, features = fe)+ RotatedAxis()+theme(legend.position = "none")+labs(title = i)
  ggsave(paste( i, "_P1.png", sep=""),width = 12,height = 10)
}

object1@reductions$umap@cell.embeddings[,1] = object1@images[["sample1"]]@coordinates[["col"]]
object1@reductions$umap@cell.embeddings[,2] = 900-object1@images[["sample1"]]@coordinates[["row"]]
plot_density(object1,features = "Plp1")+theme(axis.title = element_blank(),
                                              axis.ticks = element_blank(),
                                              axis.text.x = element_blank(),
                                              axis.text.y = element_blank(),
                                              axis.line = element_blank())

ids <- c("Erythrocytes","SLg_Fibrocyte","Out structure","Sch_Neu","Out structure","Out structure",
         "Out structure","Stria_RMC","KO",
         "Neutrophils","SLb_Fib","Corti","Macrophages")
names(ids) <- levels(object1)
object1 <- RenameIdents(object1,ids)
DimPlot(object1, reduction = "umap", label = F,pt.size=2)+NoLegend()

##细胞注释
SS8 <- subset(object1,idents=c("Sch_Neu"))
seq <- seq(0.1, 2, by = 0.1)
for(res in seq){
  sce_all <- FindClusters(SS8, resolution = seq)
}
clustree(sce_all, prefix = 'Spatial_snn_res.') + coord_flip()
SS8 <- FindClusters(SS8, verbose = FALSE,resolution =0.5)
SS8 <- RunUMAP(SS8, dims = 1:14)
DimPlot(SS8, reduction = "umap", label = T,pt.size=1)+NoLegend()

celltype_1 <- colnames(pnas_marker)
for (i in celltype_1){
  fe <- pnas_marker[i]
  p1 <- DotPlot(SS8, features = fe)+ RotatedAxis()+theme(legend.position = "none")+labs(title = i)
  ggsave(paste( i, "_P8.png", sep=""),width = 12,height = 10)
}
for (i in colnames(E16_marker)){
  fe <- E16_marker[i]
  p1 <- DotPlot(SS8, features = fe)+ RotatedAxis()+theme(legend.position = "none")+labs(title = i)
  ggsave(paste( i, "_E16.png", sep=""),width = 12,height = 10)
}


for (i in colnames(sc_top100)){
  fe <- sc_top100[i]
  p1 <- DotPlot(SS8, features = fe)+ RotatedAxis()+theme(legend.position = "none")+labs(title = i)
  ggsave(paste( i, "_adult.png", sep=""),width = 12,height = 10)
}
for (i in colnames(HS_marker)){
  fe <- HS_marker[i]
  p1 <- DotPlot(SS8, features = fe)+ RotatedAxis()+theme(legend.position = "none")+labs(title = i)
  ggsave(paste( i, "_P1.png", sep=""),width = 12,height = 10)
}

ids <- c("Neuron","Schwann cell")
names(ids) <- levels(SS8)
SS8 <- RenameIdents(SS8,ids)
Idents(object1) <- Idents(SS8)
DimPlot(object1, reduction = "umap", label = T,pt.size=1)+NoLegend()



Fib <- subset(object1,idents = c("Stria_RMC"))
seq <- seq(0.1, 2, by = 0.1)
for(res in seq){
  sce_all <- FindClusters(Fib, resolution = seq)
}
clustree(sce_all, prefix = 'Spatial_snn_res.') + coord_flip()
Fib <- FindClusters(Fib, verbose = FALSE,resolution =0.1)
Fib <- RunUMAP(Fib, dims = 1:14)
DimPlot(Fib, reduction = "umap", label = T,pt.size=1)+NoLegend()

celltype_1 <- colnames(pnas_marker)
for (i in celltype_1){
  fe <- pnas_marker[i]
  p1 <- DotPlot(Fib, features = fe)+ RotatedAxis()+theme(legend.position = "none")+labs(title = i)
  ggsave(paste( i, "_P8.png", sep=""),width = 12,height = 10)
}
for (i in colnames(E16_marker)){
  fe <- E16_marker[i]
  p1 <- DotPlot(Fib, features = fe)+ RotatedAxis()+theme(legend.position = "none")+labs(title = i)
  ggsave(paste( i, "_E16.png", sep=""),width = 12,height = 10)
}


for (i in colnames(sc_top100)){
  fe <- sc_top100[i]
  p1 <- DotPlot(Fib, features = fe)+ RotatedAxis()+theme(legend.position = "none")+labs(title = i)
  ggsave(paste( i, "_adult.png", sep=""),width = 12,height = 10)
}
for (i in colnames(HS_marker)){
  fe <- HS_marker[i]
  p1 <- DotPlot(Fib, features = fe)+ RotatedAxis()+theme(legend.position = "none")+labs(title = i)
  ggsave(paste( i, "_P1.png", sep=""),width = 12,height = 10)
}

ids <- c("Stria","Reissner's membrane")
names(ids) <- levels(Fib)
Fib <- RenameIdents(Fib,ids)
Idents(object1) <- Idents(Fib)
DimPlot(object1, reduction = "umap", label = T,pt.size=1)+NoLegend()



Fib <- subset(object1,idents = c("KO"))
seq <- seq(0.1, 2, by = 0.1)
for(res in seq){
  sce_all <- FindClusters(Fib, resolution = seq)
}
clustree(sce_all, prefix = 'Spatial_snn_res.') + coord_flip()
Fib <- FindClusters(Fib, verbose = FALSE,resolution =0.4)
Fib <- RunUMAP(Fib, dims = 1:14)
DimPlot(Fib, reduction = "umap", label = T,pt.size=1)+NoLegend()
celltype_1 <- colnames(pnas_marker)
for (i in celltype_1){
  fe <- pnas_marker[i]
  p1 <- DotPlot(Fib, features = fe)+ RotatedAxis()+theme(legend.position = "none")+labs(title = i)
  ggsave(paste( i, "_P8.png", sep=""),width = 12,height = 10)
}
for (i in colnames(E16_marker)){
  fe <- E16_marker[i]
  p1 <- DotPlot(Fib, features = fe)+ RotatedAxis()+theme(legend.position = "none")+labs(title = i)
  ggsave(paste( i, "_E16.png", sep=""),width = 12,height = 10)
}


for (i in colnames(sc_top100)){
  fe <- sc_top100[i]
  p1 <- DotPlot(Fib, features = fe)+ RotatedAxis()+theme(legend.position = "none")+labs(title = i)
  ggsave(paste( i, "_adult.png", sep=""),width = 12,height = 10)
}
for (i in colnames(HS_marker)){
  fe <- HS_marker[i]
  p1 <- DotPlot(Fib, features = fe)+ RotatedAxis()+theme(legend.position = "none")+labs(title = i)
  ggsave(paste( i, "_P1.png", sep=""),width = 12,height = 10)
}

ids <- c("Kolliker's organ","Ube2c+","Kolliker's organ")
names(ids) <- levels(Fib)
Fib <- RenameIdents(Fib,ids)
Idents(object1) <- Idents(Fib)
DimPlot(object1, reduction = "umap", label = T,pt.size=1)+NoLegend()


SS2 <- subset(object1,idents = c("Corti"))
for(res in seq){
  sce_all <- FindClusters(SS2, resolution = seq)
}
clustree(sce_all, prefix = 'Spatial_snn_res.') + coord_flip()
SS2 <- FindClusters(SS2, verbose = FALSE,resolution =0.8)
SS2 <- RunUMAP(SS2, dims = 1:14)
DimPlot(SS2, reduction = "umap", label = T,pt.size=1)+NoLegend()
celltype_1 <- colnames(pnas_marker)
for (i in celltype_1){
  fe <- pnas_marker[i]
  p1 <- DotPlot(SS2, features = fe)+ RotatedAxis()+theme(legend.position = "none")+labs(title = i)
  ggsave(paste( i, "_P8.png", sep=""),width = 12,height = 10)
}
for (i in colnames(E16_marker)){
  fe <- E16_marker[i]
  p1 <- DotPlot(SS2, features = fe)+ RotatedAxis()+theme(legend.position = "none")+labs(title = i)
  ggsave(paste( i, "_E16.png", sep=""),width = 12,height = 10)
}


for (i in colnames(sc_top100)){
  fe <- sc_top100[i]
  p1 <- DotPlot(SS2, features = fe)+ RotatedAxis()+theme(legend.position = "none")+labs(title = i)
  ggsave(paste( i, "_adult.png", sep=""),width = 12,height = 10)
}
for (i in colnames(HS_marker)){
  fe <- HS_marker[i]
  p1 <- DotPlot(SS2, features = fe)+ RotatedAxis()+theme(legend.position = "none")+labs(title = i)
  ggsave(paste( i, "_P1.png", sep=""),width = 12,height = 10)
}

ids <- c("Out sulcus cell","Deiter cell/Pillar cell","Out sulcus cell","Hair cell")
names(ids) <- levels(SS2)
SS2 <- RenameIdents(SS2,ids)
Idents(object1) <- Idents(SS2)
DimPlot(object1, reduction = "umap", label = T,pt.size=1)+NoLegend()

saveRDS(object1,"E:\\SCT\\data\\E17.5_anno\\E17_SCT_R_anno.Rds")


###################################### E17.5_L annotation #################################

##E17.5_L
celltype_1 <- colnames(pnas_marker)
for (i in celltype_1){
  fe <- pnas_marker[i]
  p1 <- DotPlot(object2, features = fe)+ RotatedAxis()+theme(legend.position = "none")+labs(title = i)
  ggsave(paste( i, ".png", sep=""),width = 12,height = 10)
}
for (i in colnames(sc_top100)){
  fe <- sc_top100[i]
  p1 <- DotPlot(object2, features = fe)+ RotatedAxis()+theme(legend.position = "none")+labs(title = i)
  ggsave(paste( i, ".png", sep=""),width = 12,height = 10)
}

SpatialDimPlot(object2, cells.highlight = CellsByIdentities(object = object1, idents = c(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15)), facet.highlight = TRUE,cols.highlight = c('#0000E3','#EDEEEF'),ncol = 4
               ,pt.size.factor = 20,alpha = c(0.5,2),stroke = 0.1,image.alpha = 0.1,crop = T)

ids <- c("SS1","SS2","Fib","SS3","SS4","SS5","SS6","SS7","SS8",
         "corti","Neuron","Out salcus","Stia")
names(ids) <- levels(object2)
object2 <- RenameIdents(object2,ids)
DimPlot(object2, reduction = "umap", label = T,pt.size=1)+NoLegend()

##细胞注释
SS8 <- subset(object2,idents=c("SS8"))
seq <- seq(0.1, 2, by = 0.1)
for(res in seq){
  sce_all <- FindClusters(SS8, resolution = seq)
}
clustree(sce_all, prefix = 'Spatial_snn_res.') + coord_flip()
SS8 <- FindClusters(SS8, verbose = FALSE,resolution =1)
SS8 <- RunUMAP(SS8, dims = 1:14)
DimPlot(SS8, reduction = "umap", label = T,pt.size=1)+NoLegend()
SpatialDimPlot(SS8, cells.highlight = CellsByIdentities(object = SS8, idents = c(0,1,2,3,4,5)), facet.highlight = TRUE,cols.highlight = c('#0000E3','#EDEEEF'),ncol = 3
               ,pt.size.factor = 20,alpha = c(0.5,2),stroke = 0.1,image.alpha = 0.1,crop = T)

celltype_1 <- colnames(pnas_marker)
for (i in celltype_1){
  fe <- pnas_marker[i]
  p1 <- DotPlot(SS8, features = fe)+ RotatedAxis()+theme(legend.position = "none")+labs(title = i)
  ggsave(paste( i, ".png", sep=""),width = 12,height = 10)
}
for (i in colnames(sc_top100)){
  fe <- sc_top100[i]
  p1 <- DotPlot(SS8, features = fe)+ RotatedAxis()+theme(legend.position = "none")+labs(title = i)
  ggsave(paste( i, ".png", sep=""),width = 12,height = 10)
}
ids <- c("Tympanic border","SS8","Reissner membrane","SS8")
names(ids) <- levels(SS8)
SS8 <- RenameIdents(SS8,ids)
Idents(object1) <- Idents(SS8)
DimPlot(object2, reduction = "umap", label = T,pt.size=1)+NoLegend()


Fib <- subset(object2,idents = c("Fib"))
seq <- seq(0.1, 2, by = 0.1)
for(res in seq){
  sce_all <- FindClusters(Fib, resolution = seq)
}
clustree(sce_all, prefix = 'Spatial_snn_res.') + coord_flip()
Fib <- FindClusters(Fib, verbose = FALSE,resolution =0.6)
Fib <- RunUMAP(Fib, dims = 1:14)
DimPlot(Fib, reduction = "umap", label = T,pt.size=1)+NoLegend()
SpatialDimPlot(Fib, cells.highlight = CellsByIdentities(object = Fib, idents = c(0,1,2)), facet.highlight = TRUE,cols.highlight = c('#0000E3','#EDEEEF'),ncol = 3
               ,pt.size.factor = 20,alpha = c(0.5,2),stroke = 0.1,image.alpha = 0.1,crop = T)
ids <- c("LW_Fib","SL_Fib","LW_Fib")
names(ids) <- levels(Fib)
Fib <- RenameIdents(Fib,ids)
Idents(object2) <- Idents(Fib)
DimPlot(object2, reduction = "umap", label = T,pt.size=1)+NoLegend()
Fib <- subset(object2,idents = c("SL_Fib"))
seq <- seq(0.1, 2, by = 0.1)
for(res in seq){
  sce_all <- FindClusters(Fib, resolution = seq)
}
clustree(sce_all, prefix = 'Spatial_snn_res.') + coord_flip()
Fib <- FindClusters(Fib, verbose = FALSE,resolution =0.5)
Fib <- RunUMAP(Fib, dims = 1:14)
DimPlot(Fib, reduction = "umap", label = T,pt.size=1)+NoLegend()
SpatialDimPlot(Fib, cells.highlight = CellsByIdentities(object = Fib, idents = c(0,1,2)), facet.highlight = TRUE,cols.highlight = c('#0000E3','#EDEEEF'),ncol = 3
               ,pt.size.factor = 20,alpha = c(0.5,2),stroke = 0.1,image.alpha = 0.1,crop = T)

celltype_1 <- colnames(pnas_marker)
for (i in celltype_1){
  fe <- pnas_marker[i]
  p1 <- DotPlot(object2, features = fe)+ RotatedAxis()+theme(legend.position = "none")+labs(title = i)
  ggsave(paste( i, ".png", sep=""),width = 12,height = 10)
}
for (i in colnames(sc_top100)){
  fe <- sc_top100[i]
  p1 <- DotPlot(object2, features = fe)+ RotatedAxis()+theme(legend.position = "none")+labs(title = i)
  ggsave(paste( i, ".png", sep=""),width = 12,height = 10)
}
SpatialDimPlot(object2, cells.highlight = CellsByIdentities(object = object1, idents = levels(object1)), facet.highlight = TRUE,cols.highlight = c('#0000E3','#EDEEEF'),ncol = 3
               ,pt.size.factor = 20,alpha = c(0.5,2),stroke = 0.1,image.alpha = 0.1,crop = T)

SS5 <- subset(object2,idents = c("SS5"))
for(res in seq){
  sce_all <- FindClusters(SS5, resolution = seq)
}
clustree(sce_all, prefix = 'Spatial_snn_res.') + coord_flip()
SS5 <- FindClusters(SS5, verbose = FALSE,resolution =1.2)
SS5 <- RunUMAP(SS5, dims = 1:14)
DimPlot(SS5, reduction = "umap", label = T,pt.size=1)+NoLegend()
for (i in celltype_1){
  fe <- pnas_marker[i]
  p1 <- DotPlot(SS5, features = fe)+ RotatedAxis()+theme(legend.position = "none")+labs(title = i)
  ggsave(paste( i, ".png", sep=""),width = 12,height = 10)
}
for (i in colnames(sc_top100)){
  fe <- sc_top100[i]
  p1 <- DotPlot(SS5, features = fe)+ RotatedAxis()+theme(legend.position = "none")+labs(title = i)
  ggsave(paste( i, ".png", sep=""),width = 12,height = 10)
}
ids <- c("Endothelial cells","Macrophages","Neutrophils","NA","Endothelial cells","Neutrophils")
names(ids) <- levels(SS5)
SS5 <- RenameIdents(SS5,ids)
Idents(object2) <- Idents(SS5)
DimPlot(object2, reduction = "umap", label = T,pt.size=1)+NoLegend()

SS <- subset(object2,idents = c("SS1","SS4","SS7"))
ids <- c("Out structure","Out structure","Out structure")
names(ids) <- levels(SS)
SS <- RenameIdents(SS,ids)
Idents(object1) <- Idents(SS)
DimPlot(object1, reduction = "umap", label = T,pt.size=1)+NoLegend()

SS2 <- subset(object2,idents = c("SS2"))
for(res in seq){
  sce_all <- FindClusters(SS2, resolution = seq)
}
clustree(sce_all, prefix = 'Spatial_snn_res.') + coord_flip()
SS2 <- FindClusters(SS2, verbose = FALSE,resolution =0.5)
SS2 <- RunUMAP(SS2, dims = 1:14)
DimPlot(SS2, reduction = "umap", label = T,pt.size=1)+NoLegend()
for (i in celltype_1){
  fe <- pnas_marker[i]
  p1 <- DotPlot(SS2, features = fe)+ RotatedAxis()+theme(legend.position = "none")+labs(title = i)
  ggsave(paste( i, ".png", sep=""),width = 12,height = 10)
}
for (i in colnames(sc_top100)){
  fe <- sc_top100[i]
  p1 <- DotPlot(SS2, features = fe)+ RotatedAxis()+theme(legend.position = "none")+labs(title = i)
  ggsave(paste( i, ".png", sep=""),width = 12,height = 10)
}
SpatialDimPlot(SS2, cells.highlight = CellsByIdentities(object = SS2, idents = c(0,1,2)), facet.highlight = TRUE,cols.highlight = c('#0000E3','#EDEEEF'),ncol = 5
               ,pt.size.factor = 28,alpha = c(0.5,2),stroke = 0.1,image.alpha = 0.1,crop = T)

ids <- c("Fibrocytes","Fibrocytes","Erythrocytes")
names(ids) <- levels(SS2)
SS2 <- RenameIdents(SS2,ids)
Idents(object1) <- Idents(SS2)
DimPlot(object1, reduction = "umap", label = T,pt.size=1)+NoLegend()


SS6 <- subset(object2,idents = c("SS6"))
for(res in seq){
  sce_all <- FindClusters(SS6, resolution = seq)
}
clustree(sce_all, prefix = 'Spatial_snn_res.') + coord_flip()
SS6 <- FindClusters(SS6, verbose = FALSE,resolution =0.7)
SS6 <- RunUMAP(SS6, dims = 1:14)
DimPlot(SS6, reduction = "umap", label = T,pt.size=1)+NoLegend()
for (i in celltype_1){
  fe <- pnas_marker[i]
  p1 <- DotPlot(SS6, features = fe)+ RotatedAxis()+theme(legend.position = "none")+labs(title = i)
  ggsave(paste( i, ".png", sep=""),width = 12,height = 10)
}
for (i in colnames(sc_top100)){
  fe <- sc_top100[i]
  p1 <- DotPlot(SS6, features = fe)+ RotatedAxis()+theme(legend.position = "none")+labs(title = i)
  ggsave(paste( i, ".png", sep=""),width = 12,height = 10)
}
SpatialDimPlot(SS6, cells.highlight = CellsByIdentities(object = SS6, idents = c(0,1,2)), facet.highlight = TRUE,cols.highlight = c('#0000E3','#EDEEEF'),ncol = 5
               ,pt.size.factor = 28,alpha = c(0.5,2),stroke = 0.1,image.alpha = 0.1,crop = T)

ids <- c("Endothelial cells","Schwan","Schwan")
names(ids) <- levels(SS6)
SS6 <- RenameIdents(SS6,ids)
Idents(object2) <- Idents(SS6)
DimPlot(object2, reduction = "umap", label = T,pt.size=1)+NoLegend()

SS3 <- subset(object2,idents = c("SS3"))
for(res in seq){
  sce_all <- FindClusters(SS3, resolution = seq)
}
clustree(sce_all, prefix = 'Spatial_snn_res.') + coord_flip()
SS3 <- FindClusters(SS3, verbose = FALSE,resolution =0.6)
SS3 <- RunUMAP(SS3, dims = 1:14)
DimPlot(SS3, reduction = "umap", label = T,pt.size=1)+NoLegend()
for (i in celltype_1){
  fe <- pnas_marker[i]
  p1 <- DotPlot(SS3, features = fe)+ RotatedAxis()+theme(legend.position = "none")+labs(title = i)
  ggsave(paste( i, ".png", sep=""),width = 12,height = 10)
}
for (i in colnames(sc_top100)){
  fe <- sc_top100[i]
  p1 <- DotPlot(SS3, features = fe)+ RotatedAxis()+theme(legend.position = "none")+labs(title = i)
  ggsave(paste( i, ".png", sep=""),width = 12,height = 10)
}
SpatialDimPlot(SS3, cells.highlight = CellsByIdentities(object = SS3, idents = c(0,1,2)), facet.highlight = TRUE,cols.highlight = c('#0000E3','#EDEEEF'),ncol = 5
               ,pt.size.factor = 28,alpha = c(0.5,2),stroke = 0.1,image.alpha = 0.1,crop = T)

ids <- c("Out structure2","Out structure2","Out structure2")
names(ids) <- levels(SS3)
SS3 <- RenameIdents(SS3,ids)
Idents(object2) <- Idents(SS3)
DimPlot(object2, reduction = "umap", label = T,pt.size=1)+NoLegend()
levels(object2)
Stria <- subset(object2,idents = c("Stia"))
for(res in seq){
  sce_all <- FindClusters(Stria, resolution = seq)
}
clustree(sce_all, prefix = 'Spatial_snn_res.') + coord_flip()
Stria <- FindClusters(Stria, verbose = FALSE,resolution =0.8)
Stria <- RunUMAP(Stria, dims = 1:14)
DimPlot(Stria, reduction = "umap", label = T,pt.size=1)+NoLegend()
for (i in celltype_1){
  fe <- pnas_marker[i]
  p1 <- DotPlot(Stria, features = fe)+ RotatedAxis()+theme(legend.position = "none")+labs(title = i)
  ggsave(paste( i, ".png", sep=""),width = 12,height = 10)
}
for (i in colnames(sc_top100)){
  fe <- sc_top100[i]
  p1 <- DotPlot(Stria, features = fe)+ RotatedAxis()+theme(legend.position = "none")+labs(title = i)
  ggsave(paste( i, ".png", sep=""),width = 12,height = 10)
}
SpatialDimPlot(Stria, cells.highlight = CellsByIdentities(object = Stria, idents = c(0,1,2)), facet.highlight = TRUE,cols.highlight = c('#0000E3','#EDEEEF'),ncol = 5
               ,pt.size.factor = 28,alpha = c(0.5,2),stroke = 0.1,image.alpha = 0.1,crop = T)

ids <- c("Marginal stria","Intermediate stria","Basal stria")
names(ids) <- levels(Stria)
Stria <- RenameIdents(Stria,ids)
Idents(object2) <- Idents(Stria)
DimPlot(object2, reduction = "umap", label = T,pt.size=1)+NoLegend()

OSC <- subset(object2,idents = c("Out salcus"))
for(res in seq){
  sce_all <- FindClusters(OSC, resolution = seq)
}
clustree(sce_all, prefix = 'Spatial_snn_res.') + coord_flip()
OSC <- FindClusters(OSC, verbose = FALSE,resolution =0.7)
OSC <- RunUMAP(OSC, dims = 1:14)
DimPlot(OSC, reduction = "umap", label = T,pt.size=1)+NoLegend()
for (i in celltype_1){
  fe <- pnas_marker[i]
  p1 <- DotPlot(OSC, features = fe)+ RotatedAxis()+theme(legend.position = "none")+labs(title = i)
  ggsave(paste( i, ".png", sep=""),width = 12,height = 10)
}
for (i in colnames(sc_top100)){
  fe <- sc_top100[i]
  p1 <- DotPlot(OSC, features = fe)+ RotatedAxis()+theme(legend.position = "none")+labs(title = i)
  ggsave(paste( i, ".png", sep=""),width = 12,height = 10)
}
SpatialDimPlot(OSC, cells.highlight = CellsByIdentities(object = OSC, idents = c(0,1,2)), facet.highlight = TRUE,cols.highlight = c('#0000E3','#EDEEEF'),ncol = 5
               ,pt.size.factor = 28,alpha = c(0.5,2),stroke = 0.1,image.alpha = 0.1,crop = T)

corti <- subset(object1,idents = c("corti"))
for(res in seq){
  sce_all <- FindClusters(corti, resolution = seq)
}
clustree(sce_all, prefix = 'Spatial_snn_res.') + coord_flip()
corti <- FindClusters(corti, verbose = FALSE,resolution =0.2)
corti <- RunUMAP(corti, dims = 1:14)
DimPlot(corti, reduction = "umap", label = T,pt.size=1)+NoLegend()
for (i in celltype_1){
  fe <- pnas_marker[i]
  p1 <- DotPlot(corti, features = fe)+ RotatedAxis()+theme(legend.position = "none")+labs(title = i)
  ggsave(paste( i, ".png", sep=""),width = 12,height = 10)
}
for (i in colnames(sc_top100)){
  fe <- sc_top100[i]
  p1 <- DotPlot(corti, features = fe)+ RotatedAxis()+theme(legend.position = "none")+labs(title = i)
  ggsave(paste( i, ".png", sep=""),width = 12,height = 10)
}
for (i in colnames(HS_marker)){
  fe <- HS_marker[i]
  p1 <- DotPlot(corti, features = fe)+ RotatedAxis()+theme(legend.position = "none")+labs(title = i)
  ggsave(paste( i, ".png", sep=""),width = 12,height = 10)
}

SpatialDimPlot(corti, cells.highlight = CellsByIdentities(object = corti, idents = c(0,1,2)), facet.highlight = TRUE,cols.highlight = c('#0000E3','#EDEEEF'),ncol = 5
               ,pt.size.factor = 28,alpha = c(0.5,2),stroke = 0.1,image.alpha = 0.1,crop = T)

ids <- c("Lateral_KO","HS","Interdental")
names(ids) <- levels(corti)
corti <- RenameIdents(corti,ids)

HS <- subset(corti,idents = c("HS"))
for(res in seq){
  sce_all <- FindClusters(HS, resolution = seq)
}
clustree(sce_all, prefix = 'Spatial_snn_res.') + coord_flip()
HS <- FindClusters(HS, verbose = FALSE,resolution =1.1)
HS <- RunUMAP(HS, dims = 1:14)
DimPlot(HS, reduction = "umap", label = T,pt.size=1)+NoLegend()
for (i in celltype_1){
  fe <- pnas_marker[i]
  p1 <- DotPlot(HS, features = fe)+ RotatedAxis()+theme(legend.position = "none")+labs(title = i)
  ggsave(paste( i, ".png", sep=""),width = 12,height = 10)
}
for (i in colnames(sc_top100)){
  fe <- sc_top100[i]
  p1 <- DotPlot(HS, features = fe)+ RotatedAxis()+theme(legend.position = "none")+labs(title = i)
  ggsave(paste( i, ".png", sep=""),width = 12,height = 10)
}
for (i in colnames(HS_marker)){
  fe <- HS_marker[i]
  p1 <- DotPlot(HS, features = fe)+ RotatedAxis()+theme(legend.position = "none")+labs(title = i)
  ggsave(paste( i, ".png", sep=""),width = 12,height = 10)
}

SpatialDimPlot(HS, cells.highlight = CellsByIdentities(object = HS, idents = c(0,1,2)), facet.highlight = TRUE,cols.highlight = c('#0000E3','#EDEEEF'),ncol = 5
               ,pt.size.factor = 28,alpha = c(0.5,2),stroke = 0.1,image.alpha = 0.1,crop = T)
FeaturePlot(HS,features = "Tbx2")+ RotatedAxis()+theme(legend.position = "none")

ids <- c("Deiter_Pillar","OHC","IHC")
names(ids) <- levels(HS)
HS <- RenameIdents(HS,ids)
Idents(corti) <- Idents(HS)

Idents(object2) <- Idents(corti)
DimPlot(object2, reduction = "umap", label = T,pt.size=1)+NoLegend()
levels(object1)
ids <- c("Deiter_Pillar","OHC","IHC",              
         "Lateral_KO","Interdental","Marginal stria",    
         "Intermediate stria","Basal stria","Out structure2"    
         ,"Endothelial cells","Schwan","Fibrocytes"        
         ,"Erythrocytes","Out structure1","Macrophages"       
         ,"Neutrophils","NA","LW_Fib"            
         ,"SL_Fib","Tympanic border","NA"               
         ,"Reissner membrane","Neuron","Out sulcus")
names(ids) <- levels(object2)
object2 <- RenameIdents(object2,ids)
DimPlot(object2, reduction = "umap", label = T,pt.size=1)+NoLegend()

object2 <- subset(object2,idents=c("Out structure2","Out structure1","Erythrocytes","Macrophages","Endothelial cells","Neutrophils",
                                   "Fibrocytes","LW_Fib","SL_Fib",
                                   "Neuron","Schwan",
                                   "Reissner membrane",
                                   "Marginal stria", "Intermediate stria","Basal stria",
                                   "Out sulcus","Tympanic border","Lateral_KO","Interdental",
                                   "OHC","IHC",              
                                   "Deiter_Pillar"))
Idents(object2) <- factor(Idents(object2), levels = c("Out structure2","Out structure1",
                                                      "Erythrocytes","Macrophages","Endothelial cells","Neutrophils",
                                                      "Fibrocytes","LW_Fib","SL_Fib",
                                                      "Neuron","Schwan",
                                                      "Reissner membrane",
                                                      "Marginal stria", "Intermediate stria","Basal stria",
                                                      "Out sulcus","Tympanic border","Lateral_KO","Interdental",
                                                      "OHC","IHC",              
                                                      "Deiter_Pillar"))
levels(object2)
object2$celltype <- Idents(object2)
saveRDS(object2,"E17_SCT_L_anno.Rds")
###################################### P8 annotation #################################

pnas_marker <- read_xlsx("PNAS_top50_gene.xlsx")
sc_top100 <- read.csv("mmtop100_cluster_ident.csv")

###################################### P8_Up annotation #################################

celltype_1 <- colnames(pnas_marker)
for (i in celltype_1){
  fe <- pnas_marker[i]
  p1 <- DotPlot(object3, features = fe)+ RotatedAxis()+theme(legend.position = "none")+labs(title = i)
  ggsave(paste( i, "_P8.png", sep=""),width = 12,height = 10)
}
for (i in colnames(E16_marker)){
  fe <- E16_marker[i]
  p1 <- DotPlot(object3, features = fe)+ RotatedAxis()+theme(legend.position = "none")+labs(title = i)
  ggsave(paste( i, "_E16.png", sep=""),width = 12,height = 10)
}


for (i in colnames(sc_top100)){
  fe <- sc_top100[i]
  p1 <- DotPlot(object3, features = fe)+ RotatedAxis()+theme(legend.position = "none")+labs(title = i)
  ggsave(paste( i, "_adult.png", sep=""),width = 12,height = 10)
}
for (i in colnames(HS_marker)){
  fe <- HS_marker[i]
  p1 <- DotPlot(object3, features = fe)+ RotatedAxis()+theme(legend.position = "none")+labs(title = i)
  ggsave(paste( i, "_P1.png", sep=""),width = 12,height = 10)
}

semla.data <- UpdateSeuratForSemla(object3)
semla.data <- LoadImages(semla.data, verbose = FALSE)
#确认坐标
MapLabels(semla.data, column_name = "seurat_clusters", ncol = 1,
          image_use = "raw",
          pt_alpha = 0.6,
          pt_size = 1,
          pt_stroke = 0.001,
          #colors = Nature_color,
          #colors = colorRampPalette(c(rep("#f7a895",5),rep("#6a3d9a",5),rep("#a9dce6",5)))(25),
          #crop_area = c(0.33, 0.34, 0.48, 0.53)
)& theme(panel.grid.major = element_line(linetype = "dashed"), axis.text = element_text())

MapLabels(semla.data, column_name = "seurat_clusters", ncol = 1,
          #image_use = "raw",
          pt_alpha = 0.6,
          pt_size = 2,
          pt_stroke = 0.001,
          #colors = Nature_color,
          #colors = colorRampPalette(c(rep("#f7a895",5),rep("#6a3d9a",5),rep("#a9dce6",5)))(25),
          crop_area = c(0.4, 0.18, 0.65, 0.45)
)&guides(fill = guide_legend(override.aes = list(size = 5)))&theme(legend.position = "right")

##分群展示
MapLabels(semla.data, column_name = "seurat_clusters",
          #image_use = "raw",
          split_labels = T,
          ncol=6,
          pt_alpha = 0.6,
          pt_size = 0.1,
          pt_stroke = 0.001,
          #colors = Nature_color,
          #colors = colorRampPalette(c(rep("#f7a895",5),rep("#6a3d9a",5),rep("#a9dce6",5)))(25),
          crop_area = c(0.4, 0.18, 0.65, 0.45)
)&NoLegend()
ids <- c("Out structure","Out structure","Neutrophils","SLg_Fibrocyte",
         "Out structure","Schwann cell","Osteoblasts","Neuron",
         "Osteoblasts","Schwann cell","IdC_Fib",
         "Root_OSC","Osteoblasts","Tympanic border cell",
         "Reissner's membrane","Corti","Stria","Out structure"
)
names(ids) <- levels(object3)
object3 <- RenameIdents(object3,ids)
DimPlot(object3, reduction = "umap", label = T,pt.size=2)+NoLegend()

SS8 <- subset(object3,idents=c("Stria"))
seq <- seq(0.1, 2, by = 0.1)
for(res in seq){
  sce_all <- FindClusters(SS8, resolution = seq)
}
clustree(sce_all, prefix = 'Spatial_snn_res.') + coord_flip()
SS8 <- FindClusters(SS8, verbose = FALSE,resolution =0.7)
SS8 <- RunUMAP(SS8, dims = 1:14)
DimPlot(SS8, reduction = "umap", label = T,pt.size=1)+NoLegend()

celltype_1 <- colnames(pnas_marker)
for (i in celltype_1){
  fe <- pnas_marker[i]
  p1 <- DotPlot(SS8, features = fe)+ RotatedAxis()+theme(legend.position = "none")+labs(title = i)
  ggsave(paste( i, "_P8.png", sep=""),width = 12,height = 10)
}
for (i in colnames(E16_marker)){
  fe <- E16_marker[i]
  p1 <- DotPlot(SS8, features = fe)+ RotatedAxis()+theme(legend.position = "none")+labs(title = i)
  ggsave(paste( i, "_E16.png", sep=""),width = 12,height = 10)
}
for (i in colnames(sc_top100)){
  fe <- sc_top100[i]
  p1 <- DotPlot(SS8, features = fe)+ RotatedAxis()+theme(legend.position = "none")+labs(title = i)
  ggsave(paste( i, "_adult.png", sep=""),width = 12,height = 10)
}
for (i in colnames(HS_marker)){
  fe <- HS_marker[i]
  p1 <- DotPlot(SS8, features = fe)+ RotatedAxis()+theme(legend.position = "none")+labs(title = i)
  ggsave(paste( i, "_P1.png", sep=""),width = 12,height = 10)
}

semla.data <- UpdateSeuratForSemla(SS8)
semla.data <- LoadImages(semla.data, verbose = FALSE)
MapLabels(semla.data, column_name = "seurat_clusters",
          #image_use = "raw",
          split_labels = T,
          ncol=3,
          pt_alpha = 0.6,
          pt_size = 0.1,
          pt_stroke = 0.001,
          #colors = Nature_color,
          #colors = colorRampPalette(c(rep("#f7a895",5),rep("#6a3d9a",5),rep("#a9dce6",5)))(25),
          crop_area = c(0.4, 0.18, 0.65, 0.45)
)&NoLegend()

ids <- c("Marginal stria","Intermediate stria","Basal stria","Marginal stria")
names(ids) <- levels(SS8)
SS8 <- RenameIdents(SS8,ids)
Idents(object3) <- Idents(SS8)
DimPlot(object3, reduction = "umap", label = T,pt.size=1)+NoLegend()

Fib <- subset(object3,idents = c("IdC_Fib"))
seq <- seq(0.1, 2, by = 0.1)
for(res in seq){
  sce_all <- FindClusters(Fib, resolution = seq)
}
clustree(sce_all, prefix = 'Spatial_snn_res.') + coord_flip()
Fib <- FindClusters(Fib, verbose = FALSE,resolution =0.4)
Fib <- RunUMAP(Fib, dims = 1:14)
DimPlot(Fib, reduction = "umap", label = T,pt.size=1)+NoLegend()
semla.data <- UpdateSeuratForSemla(Fib)
semla.data <- LoadImages(semla.data, verbose = FALSE)
MapLabels(semla.data, column_name = "seurat_clusters",
          #image_use = "raw",
          split_labels = T,
          ncol=3,
          pt_alpha = 0.6,
          pt_size = 2,
          pt_stroke = 0.001,
          #colors = Nature_color,
          #colors = colorRampPalette(c(rep("#f7a895",5),rep("#6a3d9a",5),rep("#a9dce6",5)))(25),
          crop_area = c(0.4, 0.18, 0.65, 0.45)
)&NoLegend()

celltype_1 <- colnames(pnas_marker)
for (i in celltype_1){
  fe <- pnas_marker[i]
  p1 <- DotPlot(Fib, features = fe)+ RotatedAxis()+theme(legend.position = "none")+labs(title = i)
  ggsave(paste( i, "_P8.png", sep=""),width = 12,height = 10)
}
for (i in colnames(E16_marker)){
  fe <- E16_marker[i]
  p1 <- DotPlot(Fib, features = fe)+ RotatedAxis()+theme(legend.position = "none")+labs(title = i)
  ggsave(paste( i, "_E16.png", sep=""),width = 12,height = 10)
}

for (i in colnames(sc_top100)){
  fe <- sc_top100[i]
  p1 <- DotPlot(Fib, features = fe)+ RotatedAxis()+theme(legend.position = "none")+labs(title = i)
  ggsave(paste( i, "_adult.png", sep=""),width = 12,height = 10)
}
for (i in colnames(HS_marker)){
  fe <- HS_marker[i]
  p1 <- DotPlot(Fib, features = fe)+ RotatedAxis()+theme(legend.position = "none")+labs(title = i)
  ggsave(paste( i, "_P1.png", sep=""),width = 12,height = 10)
}

ids <- c("SLb_Fibrocyte","Inner sulcus cell","Interdental cells")
names(ids) <- levels(Fib)
Fib <- RenameIdents(Fib,ids)
Idents(object3) <- Idents(Fib)
DimPlot(object3, reduction = "umap", label = T,pt.size=1)+NoLegend()

Fib <- subset(object3,idents = c("Root_OSC"))
seq <- seq(0.1, 2, by = 0.1)
for(res in seq){
  sce_all <- FindClusters(Fib, resolution = seq)
}
clustree(sce_all, prefix = 'Spatial_snn_res.') + coord_flip()
Fib <- FindClusters(Fib, verbose = FALSE,resolution =0.3)
Fib <- RunUMAP(Fib, dims = 1:14)
DimPlot(Fib, reduction = "umap", label = T,pt.size=1)+NoLegend()
semla.data <- UpdateSeuratForSemla(Fib)
semla.data <- LoadImages(semla.data, verbose = FALSE)
MapLabels(semla.data, column_name = "seurat_clusters",
          #image_use = "raw",
          split_labels = T,
          ncol=3,
          pt_alpha = 0.6,
          pt_size = 2,
          pt_stroke = 0.001,
          #colors = Nature_color,
          #colors = colorRampPalette(c(rep("#f7a895",5),rep("#6a3d9a",5),rep("#a9dce6",5)))(25),
          crop_area = c(0.4, 0.18, 0.65, 0.45)
)&NoLegend()

celltype_1 <- colnames(pnas_marker)
for (i in celltype_1){
  fe <- pnas_marker[i]
  p1 <- DotPlot(Fib, features = fe)+ RotatedAxis()+theme(legend.position = "none")+labs(title = i)
  ggsave(paste( i, "_P8.png", sep=""),width = 12,height = 10)
}
for (i in colnames(E16_marker)){
  fe <- E16_marker[i]
  p1 <- DotPlot(Fib, features = fe)+ RotatedAxis()+theme(legend.position = "none")+labs(title = i)
  ggsave(paste( i, "_E16.png", sep=""),width = 12,height = 10)
}
for (i in colnames(sc_top100)){
  fe <- sc_top100[i]
  p1 <- DotPlot(Fib, features = fe)+ RotatedAxis()+theme(legend.position = "none")+labs(title = i)
  ggsave(paste( i, "_adult.png", sep=""),width = 12,height = 10)
}
for (i in colnames(HS_marker)){
  fe <- HS_marker[i]
  p1 <- DotPlot(Fib, features = fe)+ RotatedAxis()+theme(legend.position = "none")+labs(title = i)
  ggsave(paste( i, "_P1.png", sep=""),width = 12,height = 10)
}

ids <- c("SLg_Fibrocyte","Spiral prominence","Outer sulcus cell")
names(ids) <- levels(Fib)
Fib <- RenameIdents(Fib,ids)
Idents(object3) <- Idents(Fib)
DimPlot(object3, reduction = "umap", label = T,pt.size=1)+NoLegend()


SS2 <- subset(object1,idents = c("Corti"))
for(res in seq){
  sce_all <- FindClusters(SS2, resolution = seq)
}
clustree(sce_all, prefix = 'Spatial_snn_res.') + coord_flip()
SS2 <- FindClusters(SS2, verbose = FALSE,resolution =0.9)
SS2 <- RunUMAP(SS2, dims = 1:14)
DimPlot(SS2, reduction = "umap", label = T,pt.size=1)+NoLegend()
semla.data <- UpdateSeuratForSemla(SS2)
semla.data <- LoadImages(semla.data, verbose = FALSE)
MapLabels(semla.data, column_name = "seurat_clusters",
          #image_use = "raw",
          split_labels = T,
          ncol=5,
          pt_alpha = 0.6,
          pt_size = 2,
          pt_stroke = 0.001,
          #colors = Nature_color,
          #colors = colorRampPalette(c(rep("#f7a895",5),rep("#6a3d9a",5),rep("#a9dce6",5)))(25),
          crop_area = c(0.4, 0.18, 0.65, 0.45)
)&NoLegend()


celltype_1 <- colnames(pnas_marker)
for (i in celltype_1){
  fe <- pnas_marker[i]
  p1 <- DotPlot(SS2, features = fe)+ RotatedAxis()+theme(legend.position = "none")+labs(title = i)
  ggsave(paste( i, "_P8.png", sep=""),width = 12,height = 10)
}
for (i in colnames(E16_marker)){
  fe <- E16_marker[i]
  p1 <- DotPlot(SS2, features = fe)+ RotatedAxis()+theme(legend.position = "none")+labs(title = i)
  ggsave(paste( i, "_E16.png", sep=""),width = 12,height = 10)
}


for (i in colnames(sc_top100)){
  fe <- sc_top100[i]
  p1 <- DotPlot(SS2, features = fe)+ RotatedAxis()+theme(legend.position = "none")+labs(title = i)
  ggsave(paste( i, "_adult.png", sep=""),width = 12,height = 10)
}
for (i in colnames(HS_marker)){
  fe <- HS_marker[i]
  p1 <- DotPlot(SS2, features = fe)+ RotatedAxis()+theme(legend.position = "none")+labs(title = i)
  ggsave(paste( i, "_P1.png", sep=""),width = 12,height = 10)
}

ids <- c("Deiter cell/Pillar cell","Inner Phalangeal cell","Hair cell","Inner border cell","Inner Phalangeal cell")
names(ids) <- levels(SS2)
SS2 <- RenameIdents(SS2,ids)
Idents(object3) <- Idents(SS2)
DimPlot(object3, reduction = "umap", label = T,pt.size=1)+NoLegend()
object1$celltype <- Idents(object3)
semla.data <- UpdateSeuratForSemla(object3)
semla.data <- LoadImages(semla.data, verbose = FALSE)

MapLabels(semla.data, column_name = "celltype", ncol = 1,
          #image_use = "raw",
          pt_alpha = 0.6,
          pt_size = 2,
          pt_stroke = 0.001,
          #colors = Nature_color,
          #colors = colorRampPalette(c(rep("#f7a895",5),rep("#6a3d9a",5),rep("#a9dce6",5)))(25),
          crop_area = c(0.4, 0.18, 0.65, 0.45)
)&guides(fill = guide_legend(override.aes = list(size = 5)))&theme(legend.position = "right")

saveRDS(object3,"P8_SCT_up_anno.Rds")

##################################### P8_down annotation ########################################
celltypes <- colnames(sc_top100)

for (i in celltypes){
  fe <- sc_top100[i]
  p1 <- DotPlot(object4, features = fe)+ RotatedAxis()+theme(legend.position = "none")+labs(title = i)
  ggsave(paste( i, ".png", sep=""),width = 12,height = 10)
}
object4 <- subset(object4, idents=c(1,7,10))
seq <- seq(0.1, 2, by = 0.1)
for(res in seq){
  sce_all <- FindClusters(object4, resolution = seq)
}
clustree(sce_all, prefix = 'SCT_snn_res.') + coord_flip()
object4 <- FindClusters(object4, verbose = FALSE,resolution = 2)
object4 <- RunUMAP(object4, dims = 1:30)
DimPlot(object4, reduction = "umap", label = T,pt.size=1)+NoLegend()
for (i in celltypes){
  fe <- sc_top100[i]
  p1 <- DotPlot(object4, features = fe)+ RotatedAxis()+theme(legend.position = "none")+labs(title = i)
  ggsave(paste( i, ".png", sep=""),width = 12,height = 10)
}
colnames(pnas_marker)
celltype_1 <- c("Tympanic border cells","Scala vestibuli border cells","Supporting cells","Inner hair cells",
                "Outer hair cells","Interdental cells","Deiter's cells","Pillar cells","Claudius _ Inner sulcus _ Outer suclus cells",
                "Inner border _ Inner phalangeal _ Hensen's cells","Root cells","Spindle cells","Basal stria")
for (i in celltype_1){
  fe <- pnas_marker[i]
  p1 <- DotPlot(object4, features = fe)+ RotatedAxis()+theme(legend.position = "none")+labs(title = i)
  ggsave(paste( i, ".png", sep=""),width = 12,height = 10)
}
SpatialDimPlot(object4, cells.highlight = CellsByIdentities(object = object4, idents = c(0,1,2,3,4,5,6,7)), facet.highlight = TRUE,cols.highlight = c('#0000E3','#EDEEEF'),pt.size.factor = 0.01,alpha = c(0.4,1),stroke = 0.01,image.alpha = 0.1,ncol = 4)

semla.data <- UpdateSeuratForSemla(object4)
cols <- c("white",'#C408E8','#0000E3')
MapFeatures(semla.data,features = c("Emilin2"), 
            #scale = "free",
            pt_alpha = 0.5,
            pt_size = 2,
            #override_plot_dims = TRUE,
            colors = cols,
            crop_area = c(0.4, 0.5, 0.7, 0.9)) & 
  theme(plot.title = element_blank(), legend.position = "none")
object4_sub <- subset(object4,idents=c(6,3,5,7,2))
seq <- seq(0.1, 2, by = 0.1)
for(res in seq){
  sce_all <- FindClusters(object4_sub, resolution = seq)
}
clustree(sce_all, prefix = 'SCT_snn_res.') + coord_flip()
object4_sub <- FindClusters(object4_sub, verbose = FALSE,resolution = 2)
object4_sub <- RunUMAP(object4_sub, dims = 1:30)
DimPlot(object4_sub, reduction = "umap", label = T,pt.size=1)+NoLegend()
celltype_2 <- c("Tympanic border cells","Scala vestibuli border cells","Supporting cells","Inner hair cells",
                "Outer hair cells","Interdental cells","Deiter's cells","Pillar cells","Claudius _ Inner sulcus _ Outer suclus cells",
                "Inner border _ Inner phalangeal _ Hensen's cells","Root cells","Spindle cells")

for (i in celltype_2){
  fe <- pnas_marker[i]
  p1 <- DotPlot(object4_sub, features = fe)+ RotatedAxis()+theme(legend.position = "none")+labs(title = i)
  ggsave(paste( i, ".png", sep=""),width = 12,height = 10)
}
SpatialDimPlot(object4_sub, cells.highlight = CellsByIdentities(object = P8_2_main_sub, idents = c(0,1,2,3,4,5,6,7,8,9)), facet.highlight = TRUE,cols.highlight = c('#0000E3','#EDEEEF'),pt.size.factor = 0.006,alpha = c(0.4,1),stroke = 0.01,image.alpha = 0.1,ncol = 5,crop = F)

object4_HS <- subset(object4_sub,idents=c(3,7))

seq <- seq(0.1, 2, by = 0.1)
for(res in seq){
  sce_all <- FindClusters(object4_HS, resolution = seq)
}
clustree(sce_all, prefix = 'SCT_snn_res.') + coord_flip()

object4_HS <- FindClusters(object4_HS, verbose = FALSE,resolution = 1.6)
object4_HS <- RunUMAP(object4_HS, dims = 1:30)
DimPlot(object4_HS, reduction = "umap", label = T,pt.size=1)+NoLegend()
SpatialDimPlot(object4_HS, cells.highlight = CellsByIdentities(object = object4_HS, idents = c(0,"OHC",2,3,4,5)), facet.highlight = TRUE,cols.highlight = c('#0000E3','#EDEEEF'),pt.size.factor = 0.005,alpha = c(0.4,1),stroke = 0.01,image.alpha = 0.1,ncol = 5,crop = F)

for (i in celltype_2){
  fe <- pnas_marker[i]
  p1 <- DotPlot(object4_HS, features = fe)+ RotatedAxis()+theme(legend.position = "none")+labs(title = i)
  ggsave(paste( i, ".png", sep=""),width = 12,height = 10)
}
FeaturePlot(object4_HS,features = "Clic6")+ RotatedAxis()+theme(legend.position = "none")

neu_marker <- c("Meg3","Cend1","Scg2")
sch_marker <- c("Pllp","Mag","Ncmap")
Root_marker <- c("Scg2")
Spd_marker <- c("Anxa1")
Fib_marker <-c("Crym")
Rem_marker <- c("Dkk2","Nnat","Ttr")
Tym_makrer <- c("Notum","Emilin2")
Bas_marker <-c("Fat2","Hpse","Enpep")
Int_marker <- c("Otoa")
In_sucus_marker <- c("Clic6","Tecta")
pillar_deiter <- c("Hs6st2","Sox2")
hair_marker <- c("Smpx","Ocm")
Out_hair <- subset(object4_HS,idents=c(1))
ids <- c("OHC")
names(ids) <- levels(Out_hair)
Out_hair <- RenameIdents(Out_hair, ids)
Idents(object4_HS) <- Idents(Out_hair)
Idents(object4_sub) <- Idents(Out_hair)
Iner_main <- subset(object4_HS,idents=c(2))
seq <- seq(0.1, 2, by = 0.1)
for(res in seq){
  sce_all <- FindClusters(Iner_main, resolution = seq)
}
clustree(sce_all, prefix = 'SCT_snn_res.') + coord_flip()
Iner_main <- FindClusters(Iner_main, verbose = FALSE,resolution = 1.2)
Iner_main <- RunUMAP(Iner_main, dims = 1:30)
DimPlot(Iner_main, reduction = "umap", label = T,pt.size=1)+NoLegend()
SpatialDimPlot(Iner_main, cells.highlight = CellsByIdentities(object = Iner_main, idents = c(0,1,2)), facet.highlight = TRUE,cols.highlight = c('#0000E3','#EDEEEF'),pt.size.factor = 0.005,alpha = c(0.4,1),stroke = 0.01,image.alpha = 0.1,ncol = 3,crop = F)

for (i in celltype_2){
  fe <- pnas_marker[i]
  p1 <- DotPlot(Iner_main, features = fe)+ RotatedAxis()+theme(legend.position = "none")+labs(title = i)
  ggsave(paste( i, ".png", sep=""),width = 12,height = 10)
}
FeaturePlot(Iner_main,features = "Epyc")+ RotatedAxis()+theme(legend.position = "none")

Idents(object4_sub) <- Idents(Iner_main)
ids <- c("IHC","Inner border _ Inner phalangeal","Pillar")
names(ids) <- levels(Iner_main)
Iner_main <- RenameIdents(Iner_main, ids)
Idents(object4n_sub)<- Idents(Iner_main)
Idents(object4_HS)<- Idents(Iner_main)
for (i in celltype_2){
  fe <- pnas_marker[i]
  p1 <- DotPlot(object4_sub, features = fe)+ RotatedAxis()+theme(legend.position = "none")+labs(title = i)
  ggsave(paste( i, ".png", sep=""),width = 12,height = 10)
}
SpatialDimPlot(object4_sub, cells.highlight = CellsByIdentities(object = P8_2_main_sub, idents = c(3,7)), facet.highlight = TRUE,cols.highlight = c('#0000E3','#EDEEEF'),pt.size.factor = 0.005,alpha = c(0.4,1),stroke = 0.01,image.alpha = 0.1,ncol = 2,crop = F)
ggsave("P8_7.pdf",width = 12,height = 10)

DimPlot(object4_sub, reduction = "umap", label = T,pt.size=1)+NoLegend()
FeaturePlot(object4_sub,features = "Tbx2")+ RotatedAxis()+theme(legend.position = "none")

sur <- subset(object4_HS,idents=c(3))
seq <- seq(0.1, 2, by = 0.1)
for(res in seq){
  sce_all <- FindClusters(sur, resolution = seq)
}
clustree(sce_all, prefix = 'SCT_snn_res.') + coord_flip()
sur <- FindClusters(sur, verbose = FALSE,resolution = 1.4)
sur <- RunUMAP(sur, dims = 1:30)
DimPlot(sur, reduction = "umap", label = T,pt.size=1)+NoLegend()
for (i in celltype_2){
  fe <- pnas_marker[i]
  p1 <- DotPlot(sur, features = fe)+ RotatedAxis()+theme(legend.position = "none")+labs(title = i)
  ggsave(paste( i, ".png", sep=""),width = 12,height = 10)
}
SpatialDimPlot(sur, cells.highlight = CellsByIdentities(object = sur, idents = c(0,1,2)), facet.highlight = TRUE,cols.highlight = c('#0000E3','#EDEEEF'),pt.size.factor = 0.005,alpha = c(0.4,1),stroke = 0.01,image.alpha = 0.1,ncol = 3,crop = F)

ids <- c("Inner border _ Inner phalangeal","Deiter's cells","Pillar")
names(ids) <- levels(sur)
sur <- RenameIdents(sur, ids)
Idents(object4_sub)<- Idents(sur)
Idents(object4_HS)<- Idents(sur)

Int <- subset(object4_HS,idents=c(4,5,0))
seq <- seq(0.1, 2, by = 0.1)
for(res in seq){
  sce_all <- FindClusters(Int, resolution = seq)
}
clustree(sce_all, prefix = 'SCT_snn_res.') + coord_flip()
Int <- FindClusters(Int, verbose = FALSE,resolution = 1)
Int <- RunUMAP(Int, dims = 1:30)
DimPlot(Int, reduction = "umap", label = T,pt.size=1)+NoLegend()
SpatialDimPlot(Int, cells.highlight = CellsByIdentities(object = Int, idents = c(0,4,5)), facet.highlight = TRUE,cols.highlight = c('#0000E3','#EDEEEF'),pt.size.factor = 0.005,alpha = c(0.4,1),stroke = 0.01,image.alpha = 0.1,ncol = 3,crop = F)

for (i in celltype_2){
  fe <- pnas_marker[i]
  p1 <- DotPlot(Int, features = fe)+ RotatedAxis()+theme(legend.position = "none")+labs(title = i)
  ggsave(paste( i, ".png", sep=""),width = 12,height = 10)
}

ids <- c("Inner sulcus","Inner sulcus","Inner sulcus")
names(ids) <- levels(Int)
Int <- RenameIdents(Int, ids)
Idents(object4_HS)<- Idents(Int)
Idents(object4_sub) <- Idents(Int)

DimPlot(object4_sub, reduction = "umap", label = T,pt.size=1)+NoLegend()
Int <- subset(object4_sub,idents=c(4,8))
seq <- seq(0.1, 2, by = 0.1)
for(res in seq){
  sce_all <- FindClusters(Int, resolution = seq)
}
clustree(sce_all, prefix = 'SCT_snn_res.') + coord_flip()
Int <- FindClusters(Int, verbose = FALSE,resolution = 1.5)
Int <- RunUMAP(Int, dims = 1:30)
DimPlot(Int, reduction = "umap", label = T,pt.size=1)+NoLegend()
SpatialDimPlot(Int, cells.highlight = CellsByIdentities(object = Int, idents = c(0,1,2,3)), facet.highlight = TRUE,cols.highlight = c('#0000E3','#EDEEEF'),pt.size.factor = 0.005,alpha = c(0.4,1),stroke = 0.01,image.alpha = 0.1,ncol = 3,crop = F)

for (i in celltype_2){
  fe <- pnas_marker[i]
  p1 <- DotPlot(Int, features = fe)+ RotatedAxis()+theme(legend.position = "none")+labs(title = i)
  ggsave(paste( i, ".png", sep=""),width = 12,height = 10)
}

ids <- c("NA","Interdental cells","Interdental cells","Interdental cells")
names(ids) <- levels(Int)
Int <- RenameIdents(Int, ids)
Idents(object4_sub) <- Idents(Int)

Tym <- subset(object4_sub,idents=c(1))
seq <- seq(0.1, 2, by = 0.1)
for(res in seq){
  sce_all <- FindClusters(Tym, resolution = seq)
}
clustree(sce_all, prefix = 'SCT_snn_res.') + coord_flip()
Tym <- FindClusters(Tym, verbose = FALSE,resolution = 1.2)
Tym <- RunUMAP(Tym, dims = 1:30)
DimPlot(Tym, reduction = "umap", label = T,pt.size=1)+NoLegend()
SpatialDimPlot(Tym, cells.highlight = CellsByIdentities(object = Tym, idents = c(0,1,2,3)), facet.highlight = TRUE,cols.highlight = c('#0000E3','#EDEEEF'),pt.size.factor = 0.005,alpha = c(0.4,1),stroke = 0.01,image.alpha = 0.1,ncol = 3,crop = F)

for (i in celltype_2){
  fe <- pnas_marker[i]
  p1 <- DotPlot(Tym, features = fe)+ RotatedAxis()+theme(legend.position = "none")+labs(title = i)
  ggsave(paste( i, ".png", sep=""),width = 12,height = 10)
}
tym_ids <- c("Tympanic border cells","Outer suclus cells","Tympanic border cells","Tympanic border cells")
names(tym_ids) <- levels(Tym)
Tym <- RenameIdents(Tym, tym_ids)
Idents(object4_sub) <- Idents(Tym)
DimPlot(object4_sub, reduction = "umap", label = T,pt.size=1)+NoLegend()
LW <- subset(object4_sub,idents=c(2))
seq <- seq(0.1, 2, by = 0.1)
for(res in seq){
  sce_all <- FindClusters(LW, resolution = seq)
}
clustree(sce_all, prefix = 'SCT_snn_res.') + coord_flip()
LW <- FindClusters(LW, verbose = FALSE,resolution = 1.8)
LW <- RunUMAP(LW, dims = 1:30)
DimPlot(LW, reduction = "umap", label = T,pt.size=1)+NoLegend()
SpatialDimPlot(LW, cells.highlight = CellsByIdentities(object = LW, idents = c(0,1,2,3,4,5,6)), facet.highlight = TRUE,cols.highlight = c('#0000E3','#EDEEEF'),pt.size.factor = 0.005,alpha = c(0.4,1),stroke = 0.01,image.alpha = 0.1,ncol = 3,crop = F)

for (i in celltype_2){
  fe <- pnas_marker[i]
  p1 <- DotPlot(LW, features = fe)+ RotatedAxis()+theme(legend.position = "none")+labs(title = i)
  ggsave(paste( i, ".png", sep=""),width = 12,height = 10)
}
LW_ids <- c("Outer suclus cells","Outer suclus cells","Root cell","Spindle cells","Spindle cells","Root cell","Outer suclus cells")
names(LW_ids) <- levels(LW)
LW <- RenameIdents(LW, LW_ids)
Idents(object4_sub) <- Idents(LW)
DimPlot(object4_sub, reduction = "umap", label = T,pt.size=1)+NoLegend()
SpatialDimPlot(object4_sub, cells.highlight = CellsByIdentities(object = P8_2_main_sub, idents = c(0,5,6,9)), facet.highlight = TRUE,cols.highlight = c('#0000E3','#EDEEEF'),pt.size.factor = 0.005,alpha = c(0.4,1),stroke = 0.01,image.alpha = 0.1,ncol = 3,crop = F)

celltype_3 <- colnames(pnas_marker)
for (i in celltype_3){
  fe <- pnas_marker[i]
  p1 <- DotPlot(object4_sub, features = fe)+ RotatedAxis()+theme(legend.position = "none")+labs(title = i)
  ggsave(paste( i, ".png", sep=""),width = 12,height = 10)
}
levels(object4_sub)
sub_ids <- c("Outer suclus cells","Root cell","Spindle cells","Tympanic border cells",
             "NA","Interdental cell","OHC","Inner sulcus","Inner border _ Inner phalangeal","Deiter's cells","Pillar","IHC",
             "Schwann cells","Pre-osteoblasts","SL_Fibrocytes","Macrophages")
names(sub_ids) <- levels(object4_sub)
object4_sub_rename <- RenameIdents(object4_sub, sub_ids)
DimPlot(object4_sub_rename, reduction = "umap", label = T,pt.size=1)+NoLegend()

Idents(object4) <- Idents(object4_sub_rename)

DimPlot(object4, reduction = "umap", label = T,pt.size=1)+NoLegend()
SpatialDimPlot(object4, cells.highlight = CellsByIdentities(object = P8_2_main, idents = c(0,1,4)), facet.highlight = TRUE,cols.highlight = c('#0000E3','#EDEEEF'),pt.size.factor = 0.005,alpha = c(0.4,1),stroke = 0.01,image.alpha = 0.1,ncol = 3,crop = F)
ggsave("P8_12.pdf",width = 12,height = 10)
bas <- subset(LW,idents=c(4))
seq <- seq(0.1, 2, by = 0.1)
for(res in seq){
  sce_all <- FindClusters(bas, resolution = seq)
}
clustree(sce_all, prefix = 'SCT_snn_res.') + coord_flip()
bas <- FindClusters(bas, verbose = FALSE,resolution = 1)
bas <- RunUMAP(bas, dims = 1:30)
DimPlot(bas, reduction = "umap", label = T,pt.size=1)+NoLegend()
SpatialDimPlot(bas, cells.highlight = CellsByIdentities(object = bas, idents = c(0,1,2)), facet.highlight = TRUE,cols.highlight = c('#0000E3','#EDEEEF'),pt.size.factor = 0.005,alpha = c(0.4,1),stroke = 0.01,image.alpha = 0.1,ncol = 3,crop = F)

for (i in celltype_3){
  fe <- pnas_marker[i]
  p1 <- DotPlot(bas, features = fe)+ RotatedAxis()+theme(legend.position = "none")+labs(title = i)
  ggsave(paste( i, ".png", sep=""),width = 12,height = 10)
}
bas_id <- c("basal cell","LW_fibrocyte","LW_fibrocyte")
names(bas_id) <- levels(bas)
bas <- RenameIdents(bas, bas_id)
DimPlot(bas, reduction = "umap", label = T,pt.size=1)+NoLegend()
Idents(object4) <- Idents(bas)
DimPlot(object4, reduction = "umap", label = T,pt.size=1)+NoLegend()

LW <- subset(object4,idents=c(0,1,2,3,11))
DimPlot(LW, reduction = "umap", label = T,pt.size=1)+NoLegend()
levels(LW)
id <- c("LW_fibrocyte","LW_fibrocyte","LW_fibrocyte","LW_fibrocyte","NA")
names(id) <- levels(LW)
LW <- RenameIdents(LW, id)
DimPlot(LW, reduction = "umap", label = T,pt.size=1)+NoLegend()
Idents(object4) <- Idents(LW)
DimPlot(object4, reduction = "umap", label = T,pt.size=1)+NoLegend()

seq <- seq(0.1, 2, by = 0.1)
for(res in seq){
  sce_all <- FindClusters(LW, resolution = seq)
}
clustree(sce_all, prefix = 'SCT_snn_res.') + coord_flip()
LW <- FindClusters(LW, verbose = FALSE,resolution = 1)
LW <- RunUMAP(LW, dims = 1:30)
DimPlot(LW, reduction = "umap", label = T,pt.size=1)+NoLegend()
SpatialDimPlot(LW, cells.highlight = CellsByIdentities(object = LW, idents = c(0,1,2,3,4)), facet.highlight = TRUE,cols.highlight = c('#0000E3','#EDEEEF'),pt.size.factor = 0.005,alpha = c(0.4,1),stroke = 0.01,image.alpha = 0.1,ncol = 3,crop = F)
ggsave("P8_13.pdf",width = 12,height = 10)
for (i in celltype_3){
  fe <- pnas_marker[i]
  p1 <- DotPlot(LW, features = fe)+ RotatedAxis()+theme(legend.position = "none")+labs(title = i)
  ggsave(paste( i, ".png", sep=""),width = 12,height = 10)
}

for (i in celltype_3){
  fe <- pnas_marker[i]
  p1 <- DotPlot(object4, features = fe)+ RotatedAxis()+theme(legend.position = "none")+labs(title = i)
  ggsave(paste( i, ".png", sep=""),width = 12,height = 10)
}
DimPlot(object4, reduction = "umap", label = T,pt.size=1)+NoLegend()
SpatialDimPlot(object4, cells.highlight = CellsByIdentities(object = P8_2, idents = c(9,6,2,11,4,12,13,0,5,8,3)), facet.highlight = TRUE,cols.highlight = c('#0000E3','#EDEEEF'),pt.size.factor = 0.005,alpha = c(0.4,1),stroke = 0.01,image.alpha = 0.1,ncol = 4,crop = F)

stra <- subset(object4,idents=c(13))
seq <- seq(0.1, 2, by = 0.1)
for(res in seq){
  sce_all <- FindClusters(stra, resolution = seq)
}
clustree(sce_all, prefix = 'SCT_snn_res.') + coord_flip()
stra <- FindClusters(stra, verbose = FALSE,resolution = 0.2)
stra <- RunUMAP(stra, dims = 1:30)
DimPlot(stra, reduction = "umap", label = T,pt.size=1)+NoLegend()
SpatialDimPlot(stra, cells.highlight = CellsByIdentities(object = stra, idents = c(0,1)), facet.highlight = TRUE,cols.highlight = c('#0000E3','#EDEEEF'),pt.size.factor = 0.005,alpha = c(0.4,1),stroke = 0.01,image.alpha = 0.1,ncol = 2,crop = F)

id <- c("Intermediate stria","Marginal stria")
names(id) <- levels(stra)
stra <- RenameIdents(stra, id)
Idents(object4) <- Idents(stra)
DimPlot(object4, reduction = "umap", label = T,pt.size=1)+NoLegend()

rei <- subset(object4,idents=c(12))
seq <- seq(0.1, 2, by = 0.1)
for(res in seq){
  sce_all <- FindClusters(rei, resolution = seq)
}
clustree(sce_all, prefix = 'SCT_snn_res.') + coord_flip()
rei <- FindClusters(rei, verbose = FALSE,resolution = 0.8)
rei <- RunUMAP(rei, dims = 1:30)
DimPlot(rei, reduction = "umap", label = T,pt.size=1)+NoLegend()
SpatialDimPlot(rei, cells.highlight = CellsByIdentities(object = rei, idents = c(0,1)), facet.highlight = TRUE,cols.highlight = c('#0000E3','#EDEEEF'),pt.size.factor = 0.005,alpha = c(0.4,1),stroke = 0.01,image.alpha = 0.1,ncol = 2,crop = F)

id <- c("Reissner membrane","Reissner membrane")
names(id) <- levels(rei)
rei <- RenameIdents(rei, id)
Idents(object4) <- Idents(rei)
DimPlot(object4, reduction = "umap", label = T,pt.size=1)+NoLegend()

Neu <- subset(object4,idents=c("Schwann cells"))
seq <- seq(0.1, 2, by = 0.1)
for(res in seq){
  sce_all <- FindClusters(Neu, resolution = seq)
}
clustree(sce_all, prefix = 'SCT_snn_res.') + coord_flip()
Neu <- FindClusters(Neu, verbose = FALSE,resolution = 0.5)
Neu <- RunUMAP(Neu, dims = 1:30)
DimPlot(Neu, reduction = "umap", label = T,pt.size=1)+NoLegend()
SpatialDimPlot(Neu, cells.highlight = CellsByIdentities(object = Neu, idents = c(0,1,2,3)), facet.highlight = TRUE,cols.highlight = c('#0000E3','#EDEEEF'),pt.size.factor = 0.005,alpha = c(0.4,1),stroke = 0.01,image.alpha = 0.1,ncol = 2,crop = F)

id <- c("Neuron","Schwann cells","Neuron","Neuron")
names(id) <- levels(Neu)
Neu <- RenameIdents(Neu, id)
Idents(object4) <- Idents(Neu)
DimPlot(object4, reduction = "umap", label = T,pt.size=1)+NoLegend()

Sch <- subset(object4,idents=c("Neuron"))
id <- c("Schwann cells")
names(id) <- levels(Sch)
Sch <- RenameIdents(Sch, id)
DimPlot(Sch, reduction = "umap", label = T,pt.size=1)+NoLegend()

Idents(object4) <- Idents(Sch)
DimPlot(object4, reduction = "umap", label = T,pt.size=1)+NoLegend()

IN <- subset(object4,idents=c(2,4,11))
seq <- seq(0.1, 2, by = 0.1)
for(res in seq){
  sce_all <- FindClusters(IN, resolution = seq)
}
clustree(sce_all, prefix = 'SCT_snn_res.') + coord_flip()
IN <- FindClusters(IN, verbose = FALSE,resolution = 0.5)
IN <- RunUMAP(IN, dims = 1:30)
DimPlot(IN, reduction = "umap", label = T,pt.size=1)+NoLegend()
for (i in celltype_3){
  fe <- pnas_marker[i]
  p1 <- DotPlot(IN, features = fe)+ RotatedAxis()+theme(legend.position = "none")+labs(title = i)
  ggsave(paste( i, ".png", sep=""),width = 12,height = 10)
}
levels(IN)
id <- c("Osteoblasts","Fibrocytes","Osteocytes")
names(id) <- levels(IN)
IN <- RenameIdents(IN, id)
Idents(object4) <- Idents(IN)
DimPlot(object4, reduction = "umap", label = T,pt.size=1)+NoLegend()

OUT <- subset(object4,idents=c(0,3,5,8))
seq <- seq(0.1, 2, by = 0.1)
for(res in seq){
  sce_all <- FindClusters(OUT, resolution = seq)
}
clustree(sce_all, prefix = 'SCT_snn_res.') + coord_flip()
OUT <- FindClusters(OUT, verbose = FALSE,resolution = 0.1)
OUT <- RunUMAP(OUT, dims = 1:30)
DimPlot(OUT, reduction = "umap", label = T,pt.size=1)+NoLegend()
for (i in celltype_3){
  fe <- pnas_marker[i]
  p1 <- DotPlot(OUT, features = fe)+ RotatedAxis()+theme(legend.position = "none")+labs(title = i)
  ggsave(paste( i, ".png", sep=""),width = 12,height = 10)
}
levels(IN)
id <- c("Out_structure","Neutrophils")
names(id) <- levels(OUT)
OUT <- RenameIdents(OUT, id)
Idents(object4) <- Idents(OUT)
DimPlot(object4, reduction = "umap", label = T,pt.size=1)+NoLegend()

object4_main_ano <- subset(object4,idents=c("Out_structure",
                                      "Macrophages", "Neutrophils",                 
                                      "Osteoblasts","Osteocytes","Pre-osteoblasts",
                                      "Fibrocytes","SL_Fibrocytes","LW_fibrocyte",                    
                                      "Neuron","Schwann cells",
                                      "Reissner membrane",             
                                      "basal cell","Intermediate stria","Marginal stria",                  
                                      "Root cell","Spindle cells",                            
                                      "Outer suclus cells",                           
                                      "Tympanic border cells",
                                      "Inner sulcus",  
                                      "Interdental cell",              
                                      "OHC","IHC",                 
                                      "Inner border _ Inner phalangeal","Deiter's cells","Pillar"))
Idents(object4_main_ano) <- factor(Idents(object4_main_ano), levels = c("Out_structure",
                                                                  "Macrophages", "Neutrophils",                 
                                                                  "Osteoblasts","Osteocytes","Pre-osteoblasts",
                                                                  "Fibrocytes","SL_Fibrocytes","LW_fibrocyte",                    
                                                                  "Neuron","Schwann cells",
                                                                  "Reissner membrane",             
                                                                  "basal cell","Intermediate stria","Marginal stria",                  
                                                                  "Root cell","Spindle cells",                            
                                                                  "Outer suclus cells",                           
                                                                  "Tympanic border cells",
                                                                  "Inner sulcus",  
                                                                  "Interdental cell",              
                                                                  "OHC","IHC",                 
                                                                  "Inner border _ Inner phalangeal","Deiter's cells","Pillar"))
DimPlot(object4_main_ano, reduction = "umap", label = T,pt.size=1)+NoLegend()
levels(object4_main_ano)
saveRDS(object4,"P8_SCT_Down_anno.Rds")
###################################### Adult annotation #################################

###################################### Adult Left annotation #################################



###################################### Spearman’s rank correlation calculation #################################

#Fig 1D/Fig S1-2C/Fig S1-3C/ 单细胞数据与空间单细胞数据注释结果的相关性系数计算
library(ggplot2)
library(Seurat)
library(corrplot)
library(pheatmap)
library(tidydr)
library(cowplot)
library(scCustomize)
library(Nebulosa)
library(sctransform)
library(scRNAtoolVis)
library(COSG)
library(VennDiagram)
library(ggVennDiagram)
library(ggplot2)

setwd("E:/SCT/data/相关性分析/")
##单细胞参考数据集
#E16
E16_scRNA <- readRDS("E16_scRNA.Rds")
ids <- c("Ube2c+",
         "RMC","SV",
         "IdC/ISC",
         "KO","L.KO",
         "IPhC","Hensen cell",
         "CC/OSC","LER",
         "HC",
         "HC","PsC","IPC")
names(ids) <- levels(E16_scRNA)
E16_scRNA <- RenameIdents(E16_scRNA,ids)
E16_scRNA$ident <- Idents(E16_scRNA)
DimPlot(E16_scRNA,label = TRUE,reduction = "umap")&NoLegend()
E16_scRNA$cellType <- Idents(E16_scRNA)

COSG_markers <- cosg(
  E16_scRNA,
  groups = 'all',
  assay ='RNA',
  slot='data',
  mu=1,
  n_genes_user=500)

celltype <- colnames(COSG_markers$names)
COSG_result_all <- data.frame(celltype=as.character(), markers=as.character(), COSG_score=as.numeric())

for (i in celltype) {
  COSG_result_new <- data.frame(
    celltype = rep(i, 500),
    markers = COSG_markers$names[[i]][1:500],  # 使用 [[i]] 而不是 $i
    COSG_score = COSG_markers$scores[[i]][1:500]  # 同样使用 [[i]]
  )
  COSG_result_all <- rbind(COSG_result_all, COSG_result_new)
}

E16_marker <- unique(COSG_result_all$markers)
saveRDS(E16_scRNA,"E16_scRNA_RScore.Rds")
write.csv(COSG_result_all,"E16_marker_COSG_score.csv")

##P7
P7_scRNA <- readRDS("P7_scRNA.Rds")
ids <- c("Glia",
         "M.KO","L.KO","M.L.KO",
         "IPhC",
         "HeC","CC/OSC",
         "OHC","IHC",
         "DC",
         "IPC/OPC","IPC/OPC")
names(ids) <- levels(P7_scRNA)
P7_scRNA <- RenameIdents(P7_scRNA,ids)
P7_scRNA$ident <- Idents(P7_scRNA)
DimPlot(P7_scRNA,label = TRUE,reduction = "umap")&NoLegend()
P7_scRNA$cellType <- Idents(P7_scRNA)
COSG_markers <- cosg(
  P7_scRNA,
  groups = 'all',
  assay ='RNA',
  slot='data',
  mu=1,
  n_genes_user=500)

celltype <- colnames(COSG_markers$names)
COSG_result_all <- data.frame(celltype=as.character(), markers=as.character(), COSG_score=as.numeric())

for (i in celltype) {
  COSG_result_new <- data.frame(
    celltype = rep(i, 500),
    markers = COSG_markers$names[[i]][1:500],  # 使用 [[i]] 而不是 $i
    COSG_score = COSG_markers$scores[[i]][1:500]  # 同样使用 [[i]]
  )
  COSG_result_all <- rbind(COSG_result_all, COSG_result_new)
}

P7_marker <- unique(COSG_result_all$markers)
saveRDS(P7_scRNA,"P7_scRNA_RScore.Rds")
write.csv(COSG_result_all,"P7_marker_COSG_score.csv")

#Adult
Adult_scRNA <- readRDS("Adult_scRNA.Rds")
ids <- c("PeC","MaC","EnC",
         "FC1","FC2",
         "Ia-SGNs","Ib-SGNs","Ic-SGNs","II-SGNs","III-SGNs",
         "I-SchC","II-SchC","III-SchC",
         "I-SGC","II-SGC",
         "RMC",
         "BS","MS","IS",
         "IPhC/HeC",
         "TBC","CC/OSC",
         "OHC","IHC",
         "PC")
names(ids) <- levels(Adult_scRNA)
Adult_scRNA <- RenameIdents(Adult_scRNA,ids)
Adult_scRNA$ident <- Idents(Adult_scRNA)
DimPlot(Adult_scRNA,label = TRUE,reduction = "umap")&NoLegend()
Adult_scRNA$cellType <- Idents(Adult_scRNA)
COSG_markers <- cosg(
  Adult_scRNA,
  groups = 'all',
  assay ='RNA',
  slot='data',
  mu=1,
  n_genes_user=500)

celltype <- colnames(COSG_markers$names)
COSG_result_all <- data.frame(celltype=as.character(), markers=as.character(), COSG_score=as.numeric())

for (i in celltype) {
  COSG_result_new <- data.frame(
    celltype = rep(i, 500),
    markers = COSG_markers$names[[i]][1:500],  # 使用 [[i]] 而不是 $i
    COSG_score = COSG_markers$scores[[i]][1:500]  # 同样使用 [[i]]
  )
  COSG_result_all <- rbind(COSG_result_all, COSG_result_new)
}

Adult_marker <- unique(COSG_result_all$markers)
saveRDS(Adult_scRNA,"Adult_scRNA_RScore.Rds")
write.csv(COSG_result_all,"Adult_marker_COSG_score.csv")

##空间数据集调整
object1 <- readRDS("E_L2_ano_noNA_new.Rds")

ids <- c("OS","EryC","MaC","EnC","NeC",
         "FC","SLg_FC","SLb_FC",
         "SGN",
         "SchC","SchC","SchC",
         "RMC",
         "BS","MS","IS",
         "IdC","KO",
         "CC/OSC","TBC",
         "HC","HC",
         "DC/PC")
names(ids) <- levels(object1)
object1 <- RenameIdents(object1,ids)
object1$ident <- Idents(object1)
DimPlot(object1,label = TRUE,reduction = "umap")&NoLegend()
object1$cellType <- Idents(object1)
COSG_markers <- cosg(
  object1,
  groups = 'all',
  assay ='Spatial',
  slot='data',
  mu=1,
  n_genes_user=500)

celltype <- colnames(COSG_markers$names)
COSG_result_all <- data.frame(celltype=as.character(), markers=as.character(), COSG_score=as.numeric())

for (i in celltype) {
  COSG_result_new <- data.frame(
    celltype = rep(i, 500),
    markers = COSG_markers$names[[i]][1:500],  # 使用 [[i]] 而不是 $i
    COSG_score = COSG_markers$scores[[i]][1:500]  # 同样使用 [[i]]
  )
  COSG_result_all <- rbind(COSG_result_all, COSG_result_new)
}

E17_marker <- unique(COSG_result_all$markers)
saveRDS(object1,"E17_spatial_RNA_RScore.Rds")
write.csv(COSG_result_all,"E17_spatial_marker_COSG_score.csv")

object2 <- readRDS("P8_L2_anno_main_new.Rds")
ids <- c("OS","NeC","ChC","PeC","OsC",
         "SLg_FC","SLg_FC","SLg_FC","SLg_FC","SLb_FC",
         "SGN",
         "SchC","SchC","SchC",
         "RMC",
         "BS","IS","MS","SP",
         "IdC","ISC","IPhC","HeC",
         "CC/OSC","TBC",
         "OHC","IHC",
         "DC/PC")
names(ids) <- levels(object2)
object2 <- RenameIdents(object2,ids)
object2$ident <- Idents(object2)
DimPlot(object2,label = TRUE,reduction = "umap")&NoLegend()
object2$cellType <- Idents(object2)
COSG_markers <- cosg(
  object2,
  groups = 'all',
  assay ='Spatial',
  slot='data',
  mu=1,
  n_genes_user=500)

celltype <- colnames(COSG_markers$names)
COSG_result_all <- data.frame(celltype=as.character(), markers=as.character(), COSG_score=as.numeric())

for (i in celltype) {
  COSG_result_new <- data.frame(
    celltype = rep(i, 500),
    markers = COSG_markers$names[[i]][1:500],  # 使用 [[i]] 而不是 $i
    COSG_score = COSG_markers$scores[[i]][1:500]  # 同样使用 [[i]]
  )
  COSG_result_all <- rbind(COSG_result_all, COSG_result_new)
}

P8_marker <- unique(COSG_result_all$markers)
saveRDS(object2,"P8_spatial_RNA_RScore.Rds")
write.csv(COSG_result_all,"P8_spatial_marker_COSG_score.csv")

object3 <- readRDS("E:/SCT/SCT_data/Adult_L2_ano_new.Rds")
ids <- c("OS","PeC","NeC","OsC","EnC",
         "SLg_FC","SLb_FC",
         "SGN",
         "SchC","SchC","SchC",
         "RMC",
         "BS","MS","SP",
         "IdC","HeC","TBC",
         "OHC","IHC","DC/PC")
names(ids) <- levels(object3)
object3 <- RenameIdents(object3,ids)
object3$ident <- Idents(object3)
DimPlot(object3,label = TRUE,reduction = "umap")&NoLegend()
object3$cellType <- Idents(object3)
COSG_markers <- cosg(
  object3,
  groups = 'all',
  assay ='Spatial',
  slot='data',
  mu=1,
  n_genes_user=500)

celltype <- colnames(COSG_markers$names)
COSG_result_all <- data.frame(celltype=as.character(), markers=as.character(), COSG_score=as.numeric())

for (i in celltype) {
  COSG_result_new <- data.frame(
    celltype = rep(i, 500),
    markers = COSG_markers$names[[i]][1:500],  # 使用 [[i]] 而不是 $i
    COSG_score = COSG_markers$scores[[i]][1:500]  # 同样使用 [[i]]
  )
  COSG_result_all <- rbind(COSG_result_all, COSG_result_new)
}

Adult_marker <- unique(COSG_result_all$markers)
saveRDS(object3,"Adult_spatial_RNA_RScore.Rds")
write.csv(COSG_result_all,"Adult_spatial_marker_COSG_score.csv")


##计算相关性
cortable <- function(ref=list(c(pbmc_small,"RNA","groups","label1")),
                     que=list(c(pbmc_small,"RNA","groups","label2")),
                     features=NULL,overlap="false",corr.method = 'spearman') {
  ##Part 1 相关矩阵计算
  ###分组计算平均表达量
  ###FUN1 AVERAGEEXPRESSION 
  preave <- function(obj,
                     features=c(rownames(pbmc_small)),
                     group="groups",
                     assay="RNA",
                     label="label1") {
    ob1 <- obj
    DefaultAssay(ob1) <- assay
    grp1 <- group
    Idents(ob1)<-unlist(ob1[[grp1]])
    ave1 <- AverageExpression(ob1,features = features,assays=assay)
    ave1 <- ave1[[1]]
    colnames(ave1) <- paste0(label,"_",colnames(ave1))
    
    Sp1 = ave1
    avg = rowMeans(Sp1)
    Sp1 = sweep(Sp1,1,avg,"/")
    rm(avg)
    
    Sp1[is.nan(Sp1)] <- 0
    ave1 <- Sp1
    return(ave1)
  }
  ###获取所有基因
  ###FUN2 GET FEATURES
  getfeatures <- function(ref=list(c(pbmc_small,"RNA","groups","label1")),
                          que=list(c(pbmc_small,"RNA","groups","label2"))) {
    lenref <- length(ref)
    lenque <- length(que)
    
    reflist <- list()
    for (i in 1:lenref) {
      tmp <- ref[[i]][[1]]
      DefaultAssay(tmp) <- ref[[i]][[2]]
      reflist[[i]] <- tmp
    }
    geneset1 <- lapply(reflist[1:lenref],rownames)
    gene1 <- Reduce(intersect, geneset1)
    
    quelist <- list()
    for (i in 1:lenque) {
      tmp <- que[[i]][[1]]
      DefaultAssay(tmp) <- que[[i]][[2]]
      
      quelist[[i]] <- tmp
    }
    geneset2 <- lapply(quelist[1:lenque],rownames)
    gene2 <- Reduce(intersect, geneset2)
    
    final_featues <- intersect(gene1,gene2)
    return(final_featues)
  }
  if (is.null(features)) {
    features <- getfeatures(ref,que)
  }else{
    features1 <- getfeatures(ref,que)
    features <- intersect(features,features1)
  }
  ###合并2个大组的数据
  ###FUN3 MERGE AVERAGEEXPRESSION
  preave_list <- function(inlist, features=NULL) {
    len <- length(inlist)
    matlist <- list()
    for (i in 1:length(inlist)) {
      matlist1 <- preave(obj=inlist[[i]][[1]],
                         features=features,
                         assay=inlist[[i]][[2]],
                         group=inlist[[i]][[3]],
                         label=inlist[[i]][[4]])
      matlist[[i]] <- matlist1
    }
    mat <- Reduce(cbind,matlist)
    return(mat)
    #lapply(inlist,preave(),features=features,group=)
  }
  
  ###Compute cor values
  Sp1 = preave_list(ref,features=features)
  colnames(Sp1) <- paste0(colnames(Sp1),"_ref")
  Sp2 = preave_list(que,features=features)
  colnames(Sp2) <- paste0(colnames(Sp2),"_que")
  
  geTable = merge(Sp1,Sp2, by='row.names', all=F)
  ###计算相关性系数
  rownames(geTable) = geTable$Row.names
  geTable = geTable[,2:ncol(geTable)]
  # corr.method = c('spearman', 'pearson') etc.
  corr.method <- corr.method
  Corr.Coeff.Table = cor(geTable,method=corr.method)
  
  ##Part 2 显著性检验
  ###Estimate p-value
  nPermutations = 1000
  shuffled.cor.list = list()
  pb   <- txtProgressBar(1, 100, style=3)
  
  for (i in 1:nPermutations){
    shuffled = apply(geTable[,1:ncol(Sp1)],1,sample)
    shuffled2 = apply(geTable[,(ncol(Sp1)+1):ncol(geTable)],1,sample)
    shuffled = cbind(t(shuffled),t(shuffled2))
    shuffled.cor = cor(shuffled,method=corr.method)
    shuffled.cor.list[[i]] = shuffled.cor
    rm(list=c('shuffled','shuffled2','shuffled.cor'))
    if ((i %% 100) ==0){
      setTxtProgressBar(pb, (i*100)/nPermutations)
    }
  }
  
  p.value.table = matrix(ncol=ncol(geTable), nrow = ncol(geTable))
  
  rownames(p.value.table) = colnames(geTable)
  colnames(p.value.table) = colnames(geTable)
  
  shuffled.mean.table = matrix(ncol=ncol(geTable), nrow = ncol(geTable))
  rownames(shuffled.mean.table) = colnames(geTable)
  colnames(shuffled.mean.table) = colnames(geTable)
  
  a = combn(1:ncol(geTable),2)
  for (i in 1:ncol(a)){
    cor.scores = sapply(shuffled.cor.list,"[",a[1,i],a[2,i])
    shuffled.mean.table[a[1,i],a[2,i]] = mean(cor.scores)
    shuffled.mean.table[a[2,i],a[1,i]] = mean(cor.scores)
    p.value = mean(abs(cor.scores)>=abs(Corr.Coeff.Table[a[1,i],a[2,i]]))
    p.value.table[a[1,i],a[2,i]] = p.value
    p.value.table[a[2,i],a[1,i]] = p.value
    rm(list=c('cor.scores','p.value'))
    setTxtProgressBar(pb, (i*100)/ncol(a))
  }
  if (overlap=="false") {
    M <- p.value.table
    mat <- M[,grep('ref',colnames(M))]
    mat <- mat[grep('que',rownames(M)),]
    p.value.table  <- mat
    
    M <- Corr.Coeff.Table
    mat <- M[,grep('ref',colnames(M))]
    mat <- mat[grep('que',rownames(M)),]
    Corr.Coeff.Table <- mat
  }
  
  return(list(Corr.Coeff.Table,p.value.table))
}


##
corplot <- function(cor.table=NULL,
                    pva.table=NULL,
                    cutoff=0,
                    pf=NULL,
                    col=colorRampPalette(c("darkblue", "white","darkred")),
                    order="original",
                    wid=1000,
                    hei=1000,
                    label.size=1) {
  
  cor.table[cor.table<cutoff] <- 0
  if (is.null(pf)) {
    pf <- paste0("corrplot_filter",cutoff,"_",order,".jpg")
  }
  jpeg(pf,wid,hei)
  print(
    corrplot(cor.table, order=order,tl.pos="lt", method="color", tl.col="black",cl.lim=c(min(cor.table),max(cor.table)), is.corr=F,tl.cex=label.size, sig.level=(-0.05),insig="pch", pch=19, pch.cex=0.25,pch.col="black",p.mat=pva.table, col=col(200), main= paste(pf),mar=c(3,1,5,1),cl.align.text="l")
  )
  dev.off()
  
}

cor_heatmap <- function(cor.table,pva.table=NULL,
                        col= colorRampPalette(c("darkblue", "white","darkred"))(256),
                        pf=NULL,wid=1000,hei=1000,res=1200,
                        cutoff=0,
                        scale="none"){
  
  cor.table[cor.table<cutoff] <- 0
  if (is.null(pva.table)) {
    p1<-pheatmap(cor.table,scale=scale,show_colnames = T,show_rownames = T,fontsize = 10,
                 cluster_rows = F,cluster_cols = F,border_color = "NA",color = col,res=res,filename=NA)
  }else{
    pmt <- pva.table
    ssmt <- pmt< 0.01
    pmt[ssmt] <-'**'
    smt <- pmt >0.01& pmt <0.05
    pmt[smt] <- '*'
    pmt[!ssmt&!smt]<- ''
    
    p1<-pheatmap(cor.table,scale=scale,show_colnames = T,show_rownames = T,fontsize = 10,
                 cluster_rows = F,cluster_cols = F,border_color = "NA",color = col, res=res,filename=NA,
                 display_numbers = pmt,fontsize_number = 12, number_color = "black")
    
  }
  if (is.null(pf)) {
    pf <- paste0("corheatmap_filter",cutoff,".jpg")
  }
  jpeg(pf,wid,hei)
  print(p1)
  dev.off()
  
}
corplot(cor.table=Corr.Coeff.Table,pva.table=p.value.table)
cor_heatmap(Corr.Coeff.Table,p.value.table)


scRNA <- readRDS("E16_scRNA_RScore.Rds")
spRNA <- readRDS("E17_spatial_RNA_RScore.Rds")
ref_gene <- read.csv("E16_marker_COSG_score.csv")
que_gene <- read.csv("E17_spatial_marker_COSG_score.csv")
ref <- list(c(scRNA,"RNA","cellType","scRNA"))
que <- list(c(spRNA,"Spatial","cellType","spatial"))
feature <- union(ref_gene$markers,que_gene$markers)
tmp <- cortable(ref,que,features = feature)
Corr.Coeff.Table <- tmp[[1]]
p.value.table <- tmp[[2]]
write.csv(data.frame(Corr.Coeff.Table),"E16_scRNA_E17_spatial_coeff_corration.csv")
write.csv(data.frame(p.value.table),"E16_scRNA_E17_spatial_pvalue_corration.csv")


scRNA <- readRDS("P7_scRNA_RScore.Rds")
spRNA <- readRDS("P8_spatial_RNA_RScore.Rds")
ref_gene <- read.csv("P7_marker_COSG_score.csv")
que_gene <- read.csv("P8_spatial_marker_COSG_score.csv")
ref <- list(c(scRNA,"RNA","cellType","scRNA"))
que <- list(c(spRNA,"Spatial","cellType","spatial"))
ref_top100 <- ref_gene %>% group_by(celltype) %>% top_n(100,wt=COSG_score)
que_top50 <- que_gene %>% group_by(celltype) %>% top_n(50,wt=COSG_score)

feature <- union(ref_top100$markers,que_top50$markers)
tmp <- cortable(ref,que,features = feature)
Corr.Coeff.Table <- tmp[[1]]
p.value.table <- tmp[[2]]
write.csv(data.frame(Corr.Coeff.Table),"P7_scRNA_P8_spatial_coeff_corration.csv")
write.csv(data.frame(p.value.table),"P7_scRNA_P8_spatial_pvalue_corration.csv")

scRNA <- readRDS("Adult_scRNA_RScore.Rds")
spRNA <- readRDS("Adult_spatial_RNA_RScore.Rds")
ref_gene <- read.csv("Adult_marker_COSG_score.csv")
que_gene <- read.csv("Adult_spatial_marker_COSG_score.csv")
ref <- list(c(scRNA,"RNA","cellType","scRNA"))
que <- list(c(spRNA,"Spatial","cellType","spatial"))
ref_top100 <- ref_gene %>% group_by(celltype) %>% top_n(100,wt=COSG_score)
que_top50 <- que_gene %>% group_by(celltype) %>% top_n(50,wt=COSG_score)

feature <- union(ref_top100$markers,que_top50$markers)
tmp <- cortable(ref,que,features = feature)
Corr.Coeff.Table <- tmp[[1]]
p.value.table <- tmp[[2]]
write.csv(data.frame(Corr.Coeff.Table),"Adult_scRNA_Adult_spatial_coeff_corration.csv")
write.csv(data.frame(p.value.table),"Adult_scRNA_Adult_spatial_pvalue_corration.csv")

##结果可视化

cor_matrix <- read.csv("E16_scRNA_E17_spatial_coeff_corration.csv")
rownames(cor_matrix) <- cor_matrix$X
cor_matrix<- cor_matrix[,-1]
cor_matrix <- as.matrix(cor_matrix)
p <- pheatmap(t(cor_matrix),scale="none",show_colnames = T,show_rownames = T,fontsize = 10,
              cluster_rows = F,cluster_cols = F,
              border=F,border_color = "NA",
              color = colorRampPalette(c("darkblue", "white","darkred"))(256),
              cellwidth = 20,cellheight = 16,
              #gaps_row = c(1,2,3,8,11,12),
              main = "E16_scRNA VS E17_spatial",
              fontsize_number = 12, number_color = "black")

save_pheatmap_pdf <- function(x, filename, width=7, height=7) {
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  pdf(filename, width=width, height=height)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}

save_pheatmap_pdf(p,"E16_scRNA VS E17__spatial Pearson correlation coefficient.pdf",8,6)

cor_matrix <- read.csv("P7_scRNA_P8_spatial_coeff_corration.csv")
rownames(cor_matrix) <- cor_matrix$X
cor_matrix<- cor_matrix[,-1]
cor_matrix <- as.matrix(cor_matrix)
p<-pheatmap(t(cor_matrix),scale="none",show_colnames = T,show_rownames = T,fontsize = 10,
            cluster_rows = F,cluster_cols = F,
            border=F,border_color = "NA",
            color = colorRampPalette(c("darkblue", "white","darkred"))(256),
            cellwidth = 20,cellheight = 16,
            #gaps_row = c(1,2,3,8,11,12),
            main = "P7_scRNA VS P8_spatial",
            fontsize_number = 12, number_color = "black")

save_pheatmap_pdf <- function(x, filename, width=7, height=7) {
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  pdf(filename, width=width, height=height)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}

save_pheatmap_pdf(p, "P7_scRNA VS P8_spatial Pearson correlation coefficient_top50.pdf",10,6)

cor_matrix <- read.csv("Adult_scRNA_Adult_spatial_coeff_corration.csv")
rownames(cor_matrix) <- cor_matrix$X
cor_matrix<- cor_matrix[,-1]
cor_matrix <- as.matrix(cor_matrix)
p<-pheatmap(t(cor_matrix),scale="none",show_colnames = T,show_rownames = T,fontsize = 10,
            cluster_rows = F,cluster_cols = F,
            border=F,border_color = "NA",
            color = colorRampPalette(c("darkblue", "white","darkred"))(256),
            cellwidth = 20,cellheight = 16,
            #gaps_row = c(1,2,3,8,11,12),
            main = "Adult_scRNA VS Adult_spatial",
            fontsize_number = 12, number_color = "black")

save_pheatmap_pdf <- function(x, filename, width=7, height=7) {
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  pdf(filename, width=width, height=height)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}

save_pheatmap_pdf(p, "Adult_scRNA VS Adult_spatial Pearson correlation coefficient_top50.pdf",10,8)
