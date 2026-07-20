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
data <- readRDS("E:\\SCT\\SCT_data\\E_L2_ano_noNA_new.Rds")

##For E17.5
levels(data)
ids <- c("Out structure",
         "Erythrocytes","Macrophages" ,"Endothelial cells","Neutrophils","Fibrocyte",     
         "Spiral Ligament" , "Spiral limbus" , 
         "Neuron",
         "Schwann cell" ,"Schwann cell","Schwann cell",
         "Reissner's membrane",
         "Stria Vascularis","Stria Vascularis","Stria Vascularis", 
         "Spiral limbus" ,"Sensory Epithelium" ,"Sensory Epithelium","Sensory Epithelium",  
         "Sensory Epithelium","Sensory Epithelium",
         "Sensory Epithelium")
names(ids) <- levels(data)
data <- RenameIdents(data,ids)
data$region_cluster <- Idents(data)
pltdf = data.frame(data@meta.data, x=data@images[["sample1"]]@coordinates[["col"]],
                   y=900-data@images[["sample1"]]@coordinates[["row"]])
p = ggplot(pltdf, 
           aes(x, y, color=region_cluster)) + 
  geom_point( size=2) + 
  theme_bw() +
  theme_light(base_size = 15) +
  labs(title = "") +
  theme(plot.title = element_text(hjust = 0.5))

##corti分区

p+geom_rect(aes(xmin=440,xmax=477,ymin=445,ymax=493), #n1
            fill = NA,linewidth = 1,color="red") +
  geom_rect(aes(xmin=360,xmax=415,ymin=510,ymax=545), #n2
            fill = NA,linewidth = 1,color="blue") +
  geom_rect(aes(xmin=410,xmax=460,ymin=390,ymax=445), #n3
            fill = NA,linewidth = 1,color="black")

##SGN分区
p = ggplot(pltdf, 
           aes(x, y, color=region_cluster)) + 
  geom_point( size=1) + 
  theme_bw() +
  theme_light(base_size = 15) +
  labs(title = "") +
  theme(plot.title = element_text(hjust = 0.5))

p + geom_rect(aes(xmin=402,xmax=442,ymin=461,ymax=490), #S1
              fill = NA,linewidth = 1,color="red") +
  geom_rect(aes(xmin=380,xmax=402,ymin=426,ymax=453), #S2
            fill = NA,linewidth = 1,color="blue") +
  geom_rect(aes(xmin=370,xmax=398,ymin=476,ymax=508), #S3
            fill = NA,linewidth = 1,color="black")

#数据提取
#for SE
d1 = pltdf[which(pltdf$x >= 440 & pltdf$x <= 477 & pltdf$y >= 445 & pltdf$y <= 493),]
d2 = pltdf[which(pltdf$x >= 360 & pltdf$x <= 415 & pltdf$y >= 510 & pltdf$y <= 545),]
d3 = pltdf[which(pltdf$x >= 410 & pltdf$x <= 460 & pltdf$y >= 390 & pltdf$y <= 445),]
#for Neuron
d4 = pltdf[which(pltdf$x >= 402 & pltdf$x <= 442 & pltdf$y >= 461 & pltdf$y <= 490),]
d5 = pltdf[which(pltdf$x >= 380 & pltdf$x <= 402 & pltdf$y >= 426 & pltdf$y <= 453),]
d6 = pltdf[which(pltdf$x >= 370 & pltdf$x <= 398 & pltdf$y >= 476 & pltdf$y <= 508),]


##提取区域目标细胞的spot
cltp = c("Sensory Epithelium") ###change
apex = data[,colnames(data) %in% rownames(d1)[d1$region_cluster == cltp]]
apex_se_spots <- colnames(apex)

middle = data[,colnames(data) %in% rownames(d2)[d2$region_cluster == cltp]]
middle_se_spots <- colnames(middle)

basal = data[,colnames(data) %in% rownames(d3)[d3$region_cluster == cltp]]
basal_se_spots <- colnames(basal)


data$region_annotation <- data$region_cluster
##对目标spot进行命名
# 转换为字符型
data$region_annotation <- as.character(data$region_annotation)
# 赋值
data$region_annotation[colnames(data) %in% apex_se_spots] <- "Apex_SE"
data$region_annotation[colnames(data) %in% middle_se_spots] <- "Middle_SE"
data$region_annotation[colnames(data) %in% basal_se_spots] <- "Basal_SE"

##提取区域目标细胞的spot
cltp = c("Neuron") ###change
apex = data[,colnames(data) %in% rownames(d4)[d4$region_cluster == cltp]]
apex_sgn_spots <- colnames(apex)
middle = data[,colnames(data) %in% rownames(d5)[d5$region_cluster == cltp]]
middle_sgn_spots <- colnames(middle)
basal = data[,colnames(data) %in% rownames(d6)[d6$region_cluster == cltp]]
basal_sgn_spots <- colnames(basal)

data$region_annotation[colnames(data) %in% apex_sgn_spots] <- "Apex_Neuron"
data$region_annotation[colnames(data) %in% middle_sgn_spots] <- "Middle_Neuron"
data$region_annotation[colnames(data) %in% basal_sgn_spots] <- "Basal_Neuron"
#再转回因子型
data$region_annotation <- as.factor(data$region_annotation)

pltdf = data.frame(data@meta.data, x=data@images[["sample1"]]@coordinates[["col"]],
                   y=900-data@images[["sample1"]]@coordinates[["row"]])
ggplot(pltdf, 
           aes(x, y, color=region_annotation)) + 
  geom_point( size=1) + 
  theme_bw() +
  theme_light(base_size = 15) +
  labs(title = "") +
  theme(plot.title = element_text(hjust = 0.5))


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

write.csv(SE_sig_dge.all, file = "E17_SE_sig_dge_results.csv")
write.csv(Neu_sig_dge.all, file = "E17_Neu_sig_dge_results.csv")

  
#交集差异基因集进行GO/KEGG富集分析
library(clusterProfiler)
library(org.Mm.eg.db)
##for SE
GCB_c01 = SE_sig_dge.all %>% filter(avg_log2FC > 0.25)#上调基因
GCB_c02 = SE_sig_dge.all %>% filter(avg_log2FC < -0.25)#上调基因

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

write.csv(GO_result, file = "E17_SE_Apex most_GO_enrichment_results.csv", row.names = FALSE, fileEncoding = "UTF-8")
write.csv(KEGG_result, file = "E17_SE_Apex most_KEGG_enrichment_results.csv", row.names = FALSE, fileEncoding = "UTF-8")
##basal most
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

write.csv(GO_result, file = "E17_SE_basal most_GO_enrichment_results.csv", row.names = FALSE, fileEncoding = "UTF-8")
write.csv(KEGG_result, file = "E17_SE_basal most_KEGG_enrichment_results.csv", row.names = FALSE, fileEncoding = "UTF-8")


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

write.csv(GO_result, file = "E17_Neu_Apex most_GO_enrichment_results.csv", row.names = FALSE, fileEncoding = "UTF-8")
write.csv(KEGG_result, file = "E17_Neu_Apex most_KEGG_enrichment_results.csv", row.names = FALSE, fileEncoding = "UTF-8")
##basal most
GCB_c02 = Neu_sig_dge.all %>% filter(avg_log2FC < -0.25)#上调基因
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

write.csv(GO_result, file = "E17_Neu_basal most_GO_enrichment_results.csv", row.names = FALSE, fileEncoding = "UTF-8")
write.csv(KEGG_result, file = "E17_Neu_basal most_KEGG_enrichment_results.csv", row.names = FALSE, fileEncoding = "UTF-8")
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
save(GSEA_enrichment,file = "E17_Neu_GSEA_result.RData")
saveRDS(data,"E17_region_DEG.Rds")
##统计三个时期Apex Vs Basal差异基因比例
E17_SE_sig <- read.csv("E17_SE_sig_dge_results.csv")
E17_Neu_sig <- read.csv("E17_Neu_sig_dge_results.csv")
P8_SE_sig <- read.csv("P8_SE_sig_dge_results.csv")
P8_Neu_sig <- read.csv("P8_Neu_sig_dge_results.csv")
Adult_SE_sig <- read.csv("Adult_SE_sig_dge_results.csv")
Adult_Neu_sig <- read.csv("Adult_Neu_sig_dge_results.csv")

df <- data.frame(Stage=rep(c("E17.5","P8","Adult"),2),
                 Apex_most=c(nrow(E17_SE_sig %>% filter(avg_log2FC > 0.25)),
                             nrow(P8_SE_sig %>% filter(avg_log2FC > 0.25)),
                             nrow(Adult_SE_sig %>% filter(avg_log2FC > 0.25)),
                             nrow(E17_Neu_sig %>% filter(avg_log2FC > 0.25)),
                             nrow(P8_Neu_sig %>% filter(avg_log2FC > 0.25)),
                             nrow(Adult_Neu_sig %>% filter(avg_log2FC > 0.25))),
                 Basal_most=c(nrow(E17_SE_sig %>% filter(avg_log2FC < 0.25)),
                              nrow(P8_SE_sig %>% filter(avg_log2FC < 0.25)),
                              nrow(Adult_SE_sig %>% filter(avg_log2FC < 0.25)),
                              nrow(E17_Neu_sig %>% filter(avg_log2FC < 0.25)),
                              nrow(P8_Neu_sig %>% filter(avg_log2FC < 0.25)),
                              nrow(Adult_Neu_sig %>% filter(avg_log2FC < 0.25))),
                 Region=c(rep("Sensory Epithelium",3),rep("Neuron",3)))

dfbar<-data.frame(x=c(1,2,3),
                  y=df[df$Region=="Sensory Epithelium","Apex_most"],
                  stage=df[df$Region=="Sensory Epithelium","Stage"])
dfbar1<-data.frame(x=c(1,2,3),
                   y=-df[df$Region=="Sensory Epithelium","Basal_most"],
                   stage=df[df$Region=="Sensory Epithelium","Stage"])
p1 <-ggplot() +
  # 第一个柱状图
  geom_col(data = dfbar,
           mapping = aes(x = x, y = y,fill = stage),alpha=0.5) +
  geom_text(data = dfbar,mapping = aes(x = x, y = y, label = y),  # 需要x和y映射
            position = position_nudge(y = max(dfbar$y) * 0.02),
            size = 3.5)+
  # 第二个柱状图
  geom_col(data = dfbar1,
           mapping = aes(x = x, y = y,fill = stage), alpha = 0.5)+
  geom_text(data = dfbar1,mapping = aes(x = x, y = y, label = -y),  # 需要x和y映射
            position = position_nudge(y = max(dfbar$y) * 0.02),
            size = 3.5)+
  geom_col() +
  scale_fill_manual(values = c("E17.5" = "#9400D3", 
                               "P8" = "#00CED1", 
                               "Adult" = "#FF6347"))+
  geom_hline(yintercept = 0,linewidth=1)+
  
  theme_minimal()+
  theme(
    axis.title = element_text(size = 13,
                              color = "black",
                              face = "bold"),
    axis.line.y = element_line(color = "black",
                               size = 1.2,
                               arrow = arrow(length = unit(0.3, "cm"),
                                             ends="both")
                               ),
    axis.line.x = element_blank(),
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    axis.text.y = element_blank(),
    axis.title.y = element_blank(),
    panel.grid = element_blank(),
    legend.position = "top",
    legend.direction = "vertical",
    legend.justification = c(1,0),
    legend.text = element_text(size = 15)
  )+
  annotate("text", x = -Inf, y = Inf, 
           label = "Apex_most", angle = 90,
           vjust = 2, hjust = 1,size=3) +
  annotate("text", x = -Inf, y = -Inf, 
           label = "Base_most", angle = 90,
           vjust = 2, hjust = 0,size=3) +
  coord_cartesian(clip = "off")+
  labs(title = "Sensory Epithelium Apex VS Base 差异表达基因数量时期比较")

dfbar2<-data.frame(x=c(1,2,3),
                  y=df[df$Region=="Neuron","Apex_most"],
                  stage=df[df$Region=="Neuron","Stage"])
dfbar3<-data.frame(x=c(1,2,3),
                   y=-df[df$Region=="Neuron","Basal_most"],
                   stage=df[df$Region=="Neuron","Stage"])
p2 <-ggplot() +
  # 第一个柱状图
  geom_col(data = dfbar2,
           mapping = aes(x = x, y = y,fill = stage),alpha=0.5) +
  geom_text(data = dfbar2,mapping = aes(x = x, y = y, label = y),  # 需要x和y映射
            position = position_nudge(y = max(dfbar$y) * 0.02),
            size = 3.5)+
  # 第二个柱状图
  geom_col(data = dfbar3,
           mapping = aes(x = x, y = y,fill = stage), alpha = 0.5)+
  geom_text(data = dfbar3,mapping = aes(x = x, y = y, label = -y),  # 需要x和y映射
            position = position_nudge(y = max(dfbar$y) * 0.02),
            size = 3.5)+
  geom_col() +
  scale_fill_manual(values = c("E17.5" = "#9400D3", 
                               "P8" = "#00CED1", 
                               "Adult" = "#FF6347"))+
  geom_hline(yintercept = 0,linewidth=1)+
  
  theme_minimal()+
  theme(
    axis.title = element_text(size = 13,
                              color = "black",
                              face = "bold"),
    axis.line.y = element_line(color = "black",
                               size = 1.2,
                               arrow = arrow(length = unit(0.3, "cm"),
                                             ends="both")
    ),
    axis.line.x = element_blank(),
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    axis.text.y = element_blank(),
    axis.title.y = element_blank(),
    panel.grid = element_blank(),
    legend.position = "top",
    legend.direction = "vertical",
    legend.justification = c(1,0),
    legend.text = element_text(size = 15)
  )+
  annotate("text", x = -Inf, y = Inf, 
           label = "Apex_most", angle = 90,
           vjust = 2, hjust = 1,size=3) +
  annotate("text", x = -Inf, y = -Inf, 
           label = "Base_most", angle = 90,
           vjust = 2, hjust = 0,size=3) +
  coord_cartesian(clip = "off")+
  labs(title = "Neuron Apex VS Base 差异表达基因数量时期比较")
p1/p2
ggsave("Apex VS Base 差异表达基因数量时期比较.pdf",width = 5,height = 8)
# 可视化
barplot(kegg_enrich, showCategory=20)
barplot(gene_eG, showCategory=20)
dotplot(kegg_enrich, showCategory = 20)
dotplot(kegg_enrich)
emapplot(kegg_enrich)
cnetplot(gene_eG)
P8_SE_GSEA <- data.frame(GSEA_enrichment)
