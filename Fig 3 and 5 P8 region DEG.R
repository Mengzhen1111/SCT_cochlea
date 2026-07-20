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
library(clustree)
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

##读入数据
data <- readRDS("E:\\SCT\\SCT_data\\P8_L2_anno_main_new.Rds")

##分区重命名
levels(data)
##For P8
ids <- c("Out structure",
         "Neutrophils" ,"Chondrocyte", "Pericytes","Osteoblasts",      
         "Spiral Ligament" ,"Spiral Ligament" ,"Spiral Ligament" ,"Spiral Ligament" , "Spiral limbus" , 
         "Neuron",
         "Schwann cell" ,"Schwann cell","Schwann cell",
         "Reissner's membrane",
         "Stria Vascularis","Stria Vascularis", "Stria Vascularis", 
         "Root cell/Spindle cell",
         "Spiral limbus" ,"Sensory Epithelium","Sensory Epithelium",  
         "Sensory Epithelium","Sensory Epithelium","Sensory Epithelium","Sensory Epithelium","Sensory Epithelium",
         "Sensory Epithelium")
names(ids) <- levels(data)
data <- RenameIdents(data,ids)
data$region_cluster <- Idents(data)

pltdf = data.frame(data@meta.data, x=data@images[["sample1"]]@coordinates[["col"]],
                   y=900-data@images[["sample1"]]@coordinates[["row"]])

##Corti分区坐标调整
p = ggplot(pltdf, 
           aes(x, y, color=region_cluster)) + 
  geom_point( size=1) + 
  theme_bw() +
  theme_light(base_size = 15) +
  labs(title = "") +
  theme(plot.title = element_text(hjust = 0.5))


p + geom_rect(aes(xmin=465,xmax=520,ymin=78,ymax=270), #n1
              fill = NA,linewidth = 1,color="red") +
  geom_rect(aes(xmin=525,xmax=590,ymin=80,ymax=148), #n2
            fill = NA,linewidth = 1,color="black") +
  geom_rect(aes(xmin=520,xmax=590,ymin=260,ymax=330), #n3
            fill = NA,linewidth = 1,color="blue")
p + geom_rect(aes(xmin=500,xmax=540,ymin=165,ymax=219), #S1
              fill = NA,linewidth = 1,color="red") +
  geom_rect(aes(xmin=555,xmax=625,ymin=125,ymax=198), #S2
            fill = NA,linewidth = 1,color="blue") +
  geom_rect(aes(xmin=585,xmax=620,ymin=210,ymax=269), #S3
            fill = NA,linewidth = 1,color="black")

#apex,middle,basal数据提取
d1 = pltdf[which(pltdf$x >= 465 & pltdf$x <= 520 & pltdf$y >= 78 & pltdf$y <= 270),]
d2 = pltdf[which(pltdf$x >= 525 & pltdf$x <= 590 & pltdf$y >= 80 & pltdf$y <= 148),]
d3 = pltdf[which(pltdf$x >= 520 & pltdf$x <= 590 & pltdf$y >= 260 & pltdf$y <= 330),]
d4 = pltdf[which(pltdf$x >= 500 & pltdf$x <= 540 & pltdf$y >= 165 & pltdf$y <= 219),]
d5 = pltdf[which(pltdf$x >= 555 & pltdf$x <= 625 & pltdf$y >= 125 & pltdf$y <= 198),]
d6 = pltdf[which(pltdf$x >= 585 & pltdf$x <= 620 & pltdf$y >= 210 & pltdf$y <= 269),]

##提取区域目标细胞的spot
cltp = c("Sensory Epithelium") ###change
apex = data[,colnames(data) %in% rownames(d1)[d1$region_cluster == cltp]]
apex_sgn_spots <- colnames(apex)

middle = data[,colnames(data) %in% rownames(d2)[d2$region_cluster == cltp]]
middle_sgn_spots <- colnames(middle)

basal = data[,colnames(data) %in% rownames(d3)[d3$region_cluster == cltp]]
basal_sgn_spots <- colnames(basal)


data$region_annotation <- data$region_cluster
##对目标spot进行命名
# 转换为字符型
data$region_annotation <- as.character(data$region_annotation)
# 赋值
data$region_annotation[colnames(data) %in% apex_sgn_spots] <- "Apex_SE"
data$region_annotation[colnames(data) %in% middle_sgn_spots] <- "Middle_SE"
data$region_annotation[colnames(data) %in% basal_sgn_spots] <- "Basal_SE"

#Neu apex,middle,basal数据提取
cltp = c("Neuron") ###change
apex = data[,colnames(data) %in% rownames(d4)[d4$celltype == cltp]]
apex_sgn_spots <- colnames(apex)
middle = data[,colnames(data) %in% rownames(d5)[d5$celltype == cltp]]
middle_sgn_spots <- colnames(middle)
basal = data[,colnames(data) %in% rownames(d6)[d6$celltype == cltp]]
basal_sgn_spots <- colnames(basal)

data$region_annotation[colnames(data) %in% apex_sgn_spots] <- "Apex_Neuron"
data$region_annotation[colnames(data) %in% middle_sgn_spots] <- "Middle_Neuron"
data$region_annotation[colnames(data) %in% basal_sgn_spots] <- "Basal_Neuron"
#再转回因子型
data$region_annotation <- as.factor(data$region_annotation)

pltdf = data.frame(data@meta.data, x=data@images[["sample1"]]@coordinates[["col"]],
                   y=900-data@images[["sample1"]]@coordinates[["row"]])
p = ggplot(pltdf, 
           aes(x, y, color=region_annotation)) + 
  geom_point( size=1) + 
  #scale_color_manual(values = A_color)+
  theme_bw() +
  theme_light(base_size = 15) +
  labs(title = "") +
  theme(plot.title = element_text(hjust = 0.5))

p + geom_rect(aes(xmin=465,xmax=520,ymin=78,ymax=270), #n1
              fill = NA,linewidth = 1,color="red") +
  geom_rect(aes(xmin=525,xmax=590,ymin=80,ymax=148), #n2
            fill = NA,linewidth = 1,color="black") +
  geom_rect(aes(xmin=520,xmax=590,ymin=260,ymax=330), #n3
            fill = NA,linewidth = 1,color="blue")

##差异分析
#findallmarkers得到各区域特意表达的差异基因集
#Idents(data) <- data$region_cluster
#markers <- FindAllMarkers(data, only.pos = T, min.pct = 0.1, logfc.threshold = 0.1)
#write.csv(markers,"E_region_markers.csv")

#findmarkers得到特定区域比较的差异基因集
Idents(data) <- data$region_annotation
SE_AB_markers <- FindMarkers(data,
                       ident.1 ="Apex_SE",
                       ident.2="Basal_SE",assay = 'Spatial',slot = 'counts',
                       logfc.threshold =0,min.pct = 0 )

SE_sig_dge.all <- subset(SE_AB_markers, p_val < 0.05)

Neu_AB_markers <- FindMarkers(data,
                             ident.1 ="Apex_Neuron",
                             ident.2="Basal_Neuron",assay = 'Spatial',slot = 'counts',
                             logfc.threshold =0,min.pct = 0 )
Neu_sig_dge.all <- subset(Neu_AB_markers, p_val < 0.05)

write.csv(SE_sig_dge.all, file = "P8_SE_sig_dge_results.csv")
write.csv(Neu_sig_dge.all, file = "P8_Neu_sig_dge_results.csv")

  
#差异基因集进行GO/KEGG富集分析
library(clusterProfiler)
library(org.Mm.eg.db)
##for SE
GCB_c01 = SE_sig_dge.all %>% filter(avg_log2FC > 0.25)#apex most基因
gid <- bitr(unique(rownames(GCB_c01)), "SYMBOL", "ENTREZID", OrgDb = "org.Mm.eg.db")
kegg_enrich <- enrichKEGG(gene= gid$ENTREZID,keyType = 'kegg',organism = "mmu", # hsa=人类，mmu=小鼠等
                          pAdjustMethod = 'fdr',
                          pvalueCutoff= 0.1,
                          qvalueCutoff= 0.3)
kegg_enrich_symbol <- setReadable(kegg_enrich, 
                                  OrgDb= org.Mm.eg.db,
                                  keyType= "ENTREZID")
gene_eG <- enrichGO(gene = gid$ENTREZID, #需要分析的基因的EntrezID
                    OrgDb = org.Mm.eg.db, #人基因数据库
                    pvalueCutoff =0.05, #设置pvalue界值
                    qvalueCutoff = 0.2, #设置qvalue界值(FDR校正后的p值）
                    ont="all", #选择功能富集的类型，可选BP、MF、CC，这里选择all。
                    readable =T)

KEGG_result <- kegg_enrich_symbol@result
GO_result <- gene_eG@result

write.csv(GO_result, file = "P8_SE_Apex most_GO_enrichment_results.csv", row.names = FALSE, fileEncoding = "UTF-8")
write.csv(KEGG_result, file = "P8_SE_Apex most_KEGG_enrichment_results.csv", row.names = FALSE, fileEncoding = "UTF-8")
##basal most
GCB_c02 = SE_sig_dge.all %>% filter(avg_log2FC < -0.25)#basal most基因
gid <- bitr(unique(rownames(GCB_c02)), "SYMBOL", "ENTREZID", OrgDb = "org.Mm.eg.db")
kegg_enrich <- enrichKEGG(gene= gid$ENTREZID,keyType = 'kegg',organism = "mmu", # hsa=人类，mmu=小鼠等
                          pAdjustMethod = 'fdr',
                          pvalueCutoff= 0.1,
                          qvalueCutoff= 0.3)
kegg_enrich_symbol <- setReadable(kegg_enrich, 
                                  OrgDb= org.Mm.eg.db,
                                  keyType= "ENTREZID")
gene_eG <- enrichGO(gene = gid$ENTREZID, #需要分析的基因的EntrezID
                    OrgDb = org.Mm.eg.db, #人基因数据库
                    pvalueCutoff =0.05, #设置pvalue界值
                    qvalueCutoff = 0.2, #设置qvalue界值(FDR校正后的p值）
                    ont="all", #选择功能富集的类型，可选BP、MF、CC，这里选择all。
                    readable =T)

KEGG_result <- kegg_enrich_symbol@result
GO_result <- gene_eG@result

write.csv(GO_result, file = "P8_SE_basal most_GO_enrichment_results.csv", row.names = FALSE, fileEncoding = "UTF-8")
write.csv(KEGG_result, file = "P8_SE_basal most_KEGG_enrichment_results.csv", row.names = FALSE, fileEncoding = "UTF-8")



##For Neuron
GCB_c01 = Neu_sig_dge.all %>% filter(avg_log2FC > 0.25)#上调基因
gid <- bitr(unique(rownames(GCB_c01)), "SYMBOL", "ENTREZID", OrgDb = "org.Mm.eg.db")
kegg_enrich <- enrichKEGG(gene= gid$ENTREZID,keyType = 'kegg',organism = "mmu", # hsa=人类，mmu=小鼠等
                          pAdjustMethod = 'fdr',
                          pvalueCutoff= 0.1,
                          qvalueCutoff= 0.3)
kegg_enrich_symbol <- setReadable(kegg_enrich, 
                                  OrgDb= org.Mm.eg.db,
                                  keyType= "ENTREZID")
gene_eG <- enrichGO(gene = gid$ENTREZID, #需要分析的基因的EntrezID
                    OrgDb = org.Mm.eg.db, #人基因数据库
                    pvalueCutoff =0.05, #设置pvalue界值
                    qvalueCutoff = 0.2, #设置qvalue界值(FDR校正后的p值）
                    ont="all", #选择功能富集的类型，可选BP、MF、CC，这里选择all。
                    readable =T)

KEGG_result <- kegg_enrich_symbol@result
GO_result <- gene_eG@result

write.csv(GO_result, file = "P8_Neu_Apex most_GO_enrichment_results.csv", row.names = FALSE, fileEncoding = "UTF-8")
write.csv(KEGG_result, file = "P8_Neu_Apex most_KEGG_enrichment_results.csv", row.names = FALSE, fileEncoding = "UTF-8")
##basal most
GCB_c02 = Neu_sig_dge.all %>% filter(avg_log2FC < -0.25)#basal most基因
gid <- bitr(unique(rownames(GCB_c02)), "SYMBOL", "ENTREZID", OrgDb = "org.Mm.eg.db")
kegg_enrich <- enrichKEGG(gene= gid$ENTREZID,keyType = 'kegg',organism = "mmu", # hsa=人类，mmu=小鼠等
                          pAdjustMethod = 'fdr',
                          pvalueCutoff= 0.1,
                          qvalueCutoff= 0.3)
kegg_enrich_symbol <- setReadable(kegg_enrich, 
                                  OrgDb= org.Mm.eg.db,
                                  keyType= "ENTREZID")
gene_eG <- enrichGO(gene = gid$ENTREZID, #需要分析的基因的EntrezID
                    OrgDb = org.Mm.eg.db, #人基因数据库
                    pvalueCutoff =0.05, #设置pvalue界值
                    qvalueCutoff = 0.2, #设置qvalue界值(FDR校正后的p值）
                    ont="all", #选择功能富集的类型，可选BP、MF、CC，这里选择all。
                    readable =T)

KEGG_result <- kegg_enrich_symbol@result
GO_result <- gene_eG@result

write.csv(GO_result, file = "P8_Neu_basal most_GO_enrichment_results.csv", row.names = FALSE, fileEncoding = "UTF-8")
write.csv(KEGG_result, file = "P8_Neu_basal most_KEGG_enrichment_results.csv", row.names = FALSE, fileEncoding = "UTF-8")
##使用GSEA进行整体的功能富集
genelist <- Neu_AB_markers$avg_log2FC
names(genelist) <- rownames(Neu_AB_markers)
genelist <- sort(genelist, decreasing = T)
head(genelist)
GSEA_enrichment <- GSEA(genelist,                 # 排序后的gene
                        TERM2GENE = geneSet_onco, # 基因集
                        pvalueCutoff = 0.05,      # P值阈值
                        minGSSize = 20,           # 最小基因数量
                        maxGSSize = 1000,         # 最大基因数量
                        eps = 0,                  # P值边界
                        pAdjustMethod = "BH")     # 校正P值的计算方法
result <- data.frame(GSEA_enrichment)
dim(GSEA_enrichment@result)
save(GSEA_enrichment,file = "P8_Neu_GSEA_result.RData")
saveRDS(data,"P8_region_DEG.Rds")
