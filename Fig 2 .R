library(monocle)
library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(tidyverse)
library(reshape2)
library(cluster)
library(tibble)
library(celldex)
library(clustree)
library(patchwork)
library(RColorBrewer)
library(semla)
library(ggforce)
library(ggplot2)
library(cowplot)
library(scRNAtoolVis)
library(COSG)
library(readxl)
library(tidyverse)
library(RColorBrewer)
library(ComplexHeatmap)
##载入数据
object1 <- readRDS("E17_SCT_L_anno.Rds")
object2 <- readRDS("P8_SCT_Down_anno.Rds")
object3 <- readRDS("Adult_SCT_Right_anno.Rds")
object1$group <- "E17.5"
object2$group <- "P8"
object3$group <- "Adult"
##20250923
E17 <- subset(object1,idents="SGN")
E17$group <- "E17.5"
P8 <- subset(object2,idents=c("SGN"))
P8$group <- "P8"
Adult <- subset(object3,idents=c("SGN"))
Adult$group <- "Adult"
combined <- merge(E17, y = c(P8, Adult), add.cell.ids = c("E17", "P8", "Adult"))
combined <- NormalizeData(combined)
combined <- FindVariableFeatures(combined)
combined <- ScaleData(combined)
combined <- RunPCA(combined)
count1_matrix <- as.matrix(GetAssayData(combined, layer = "counts.1"))
count2_matrix <- as.matrix(GetAssayData(combined, layer = "counts.2"))
count3_matrix <- as.matrix(GetAssayData(combined, layer = "counts.3"))
genes1 <- rownames(count1_matrix)
genes2 <- rownames(count2_matrix)
genes3 <- rownames(count3_matrix)
common_genes <- Reduce(intersect, list(genes1, genes2, genes3))
count1_matrix <- count1_matrix[common_genes, ]
count2_matrix <- count2_matrix[common_genes, ]
count3_matrix <- count3_matrix[common_genes, ]
expression_matrix <- cbind(count1_matrix, count2_matrix, count3_matrix)
cell_metadata <- combined@meta.data
gene_annotation <- data.frame(gene_short_name = rownames(expression_matrix), row.names = rownames(expression_matrix))
pd <- new('AnnotatedDataFrame', data = cell_metadata)
fd <- new('AnnotatedDataFrame', data = gene_annotation)
cds <- newCellDataSet(expression_matrix, phenoData = pd, featureData = fd, 
                      expressionFamily = negbinomial.size())


##归一化处理
cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)
##数据过滤
cds<-detectGenes(cds,min_expr = 0.1)
##选择过程基因

##选择monocle计算的差异基因
diff_test_res <- differentialGeneTest(cds, fullModelFormulaStr = "~group",cores = 4)
ordering_genes <- rownames(subset(diff_test_res, pval < 0.01))

cds <- setOrderingFilter(cds,ordering_genes)
plot_ordering_genes(cds)

cds <- reduceDimension(cds,method='DDRTree')
##排序
cds <- orderCells(cds,reverse = T)
# 可视化轨迹，按伪时间着色
plot_cell_trajectory(cds, color_by = "Pseudotime")+
  tidydr::theme_dr(xlength = 0.2, ylength = 0.2,arrow = arrow(length = unit(0.2, "inches"),type = "closed")) +
  theme(panel.grid = element_blank(),
        axis.title = element_text(face = 2,hjust = 0.03))

ggsave("SCT_Neu_pseu_Pseudotime.pdf",width = 5,height = 3)

# 可视化轨迹，按样本着色
plot_cell_trajectory(cds, color_by = "group")+scale_colour_manual(values = c("#79B494", "#D67E56",  "#848CBD"))+
  tidydr::theme_dr(xlength = 0.2, ylength = 0.2,arrow = arrow(length = unit(0.2, "inches"),type = "closed")) +
  theme(panel.grid = element_blank(),
        axis.title = element_text(face = 2,hjust = 0.03))
ggsave("SCT_Neu_pseu_group.pdf",width = 5,height = 3)

plot_cell_trajectory(cds, color_by = "State")+
  tidydr::theme_dr(xlength = 0.2, ylength = 0.2,arrow = arrow(length = unit(0.2, "inches"),type = "closed")) +
  theme(panel.grid = element_blank(),
        axis.title = element_text(face = 2,hjust = 0.03))
ggsave("SCT_Neu_pseu_State.pdf",width = 5,height = 3)

cds_DGT_pseudotimegenes <- differentialGeneTest(cds, fullModelFormulaStr = "~group")
cds_DGT_pseudotimegenes_sig <- subset(cds_DGT_pseudotimegenes, qval < 0.01)

cds_DGT_pseudotimegenes_sig_noMT <- cds_DGT_pseudotimegenes_sig[!grepl("^mt-", rownames(cds_DGT_pseudotimegenes_sig)), ]

cds_DGT_pseudotimegenes_sig_noMT_Neu_marker <- intersect(cds_DGT_pseudotimegenes_sig_noMT$gene_short_name, Neu_marker)
cds_DGT_pseudotimegenes_sig_noMT_Neu_marker <- cds_DGT_pseudotimegenes_sig_noMT[cds_DGT_pseudotimegenes_sig_noMT$gene_short_name %in% cds_DGT_pseudotimegenes_sig_noMT_Neu_marker,]
write.csv(cds_DGT_pseudotimegenes_sig_noMT,"SCT_Neu_peudotimegenes_qvalue_0.01.csv")
write.csv(cds_DGT_pseudotimegenes_sig_noMT_Neu_marker,"SCT_Neu_peudotimegenes_qvalue_0.01_neuMarkers.csv")
pdf("SCT_Neu_pseudotime_336_heatmap.pdf",height = 8,width = 7)
plot_pseudotime_heatmap(cds[cds_DGT_pseudotimegenes_sig_noMT_Neu_marker$gene_short_name,], 
                        num_cluster = 8, 
                        show_rownames = F, 
                        return_heatmap = T,
                        #add_annotation_col = annotation_col,
                        hmcols = colorRampPalette(c("navy","white","firebrick3"))(100))

dev.off()
dev.new()

delibrary(clusterProfiler)
library(msigdbr)
library(GSVA)
library(org.Mm.eg.db)
library(enrichR)
View(msigdbr_collections())
setwd("../")
## 获取各分支基因

p=plot_pseudotime_heatmap(cds[cds_DGT_pseudotimegenes_sig_noMT_Neu_marker$gene_short_name,], 
                          num_cluster = 8, 
                          show_rownames = T, 
                          return_heatmap = T,
                          hmcols = colorRampPalette(c("navy","white","firebrick3"))(100))

clusters <- cutree(p$tree_row, k = 8)
clustering <- data.frame(clusters)
clustering[,1] <- as.character(clustering[,1])
colnames(clustering) <- "Gene_Clusters"
clustering$gene <- rownames(clustering)
write.csv(clustering,"SCT_Neu_peudotimegenes_qvalue_0.01_neuMarkers_subcluster.csv")
##进行趋势分类后加载数据进行GO分析
clustering <- read.csv("SCT_Neu_peudotimegenes_qvalue_0.01_neuMarkers_subcluster.csv")
##分别进行GO分析

GO_all<-data.frame(Description=character(),
                   qvalue=numeric(),
                   pvalue=numeric(),
                   RichFactor=numeric(),
                   Count=numeric(),
                   Ontology=character(),
                   gene=character(),
                   group=character())
for (i in unique(clustering$X.1)) {
  gene=clustering[clustering$X.1==i,]$gene
  gene_entrez_id2GO <- clusterProfiler::bitr(gene, "SYMBOL", "ENTREZID", "org.Mm.eg.db", drop = TRUE)$ENTREZID    
  eG <- enrichGO(gene = gene_entrez_id2GO, #需要分析的基因的EntrezID
                 OrgDb = org.Mm.eg.db, #人基因数据库
                 pvalueCutoff =0.01, #设置pvalue界值
                 qvalueCutoff = 0.01, #设置qvalue界值(FDR校正后的p值）
                 ont="all", #选择功能富集的类型，可选BP、MF、CC，这里选择all。
                 readable =T)
  GO <- data.frame(eG@result)
  GO$group <- rep(i,nrow(GO))
  GO_new <- data.frame(Description=GO$Description,
                       qvalue=GO$qvalue,
                       pvalue=GO$pvalue,
                       RichFactor=GO$RichFactor,
                       Count=GO$Count,
                       Ontology=GO$ONTOLOGY,
                       gene=GO$geneID,
                       group=GO$group)
  GO_all <- rbind(GO_all,GO_new)
}
write.csv(GO_all,"SCT_Neu_pseu_sig_NeuMarker_group_pvalue0.01_GOresult.csv")
saveRDS(cds,"SCT_Neu_pseu_new2.Rds")

clustering <- read.csv("SCT_Neu_peudotimegenes_qvalue_0.01_neuMarkers_subcluster.csv")

cluster_Adult_most <- clustering[clustering$X.1=="Adult most",]

cluster_E_most <- clustering[clustering$X.1=="E17.5 most",]
cluster_gene <- union(cluster_Adult_most$gene,cluster_E_most$gene)
plot_pseudotime_heatmap(cds[cluster_gene,], 
                        num_cluster = 2, 
                        show_rownames = F, 
                        return_heatmap = T,
                        #add_annotation_col = annotation_col,
                        hmcols = colorRampPalette(c("navy","white","firebrick3"))(100))
GO_terms <- read_xlsx("SCT_Neu_pseu_GO_displayterms.xlsx")

ggplot(GO_terms)+
  geom_bar(aes(Count,fct_reorder(Description,group),fill = Ontology),stat="identity") + #将Term按照Ontology排序
  theme_bw()+
  theme(text=element_text(family='serif',size =10,face = 'bold'),axis.text = element_text(family='serif',size =10,face = 'bold',colour = 'black')) + 
  geom_text(data=subset(GO_terms,qvalue<0.05,c('Count','Description')),aes(x=Count+3,y=as.factor(Description),label='*')) + #对于qvalue<0.05的term添加*
  labs(x='Number of genes',y='GO term') + 
  scale_fill_manual(values=c(BP = "#79B494", CC = "#D67E56", MF = "#848CBD"))

mytheme <- theme(axis.title=element_text(face="bold", size=14,colour = 'black'), #坐标轴标题
                 axis.text=element_text(face="bold", size=14,colour = 'black'), #坐标轴标签
                 axis.line = element_line(size=0.5, colour = 'black'), #轴线
                 panel.background = element_rect(color='black'), #绘图区边框
                 legend.key = element_blank() #关闭图例边框
)

ggplot(GO_terms,aes(x=as.numeric(Count),y=Description,colour=as.numeric(RichFactor),size=-1*log10(as.numeric(pvalue))))+
  geom_point()+
  scale_size(range=c(2, 8))+
  scale_colour_gradient(low = "blue",high = "red")+
  facet_grid(group ~ ., scales = "free_y", space = "free_y")+
  theme_bw()+
  ylab("GO Terms")+
  xlab("Number of gene")+
  labs(color="RichFactor",size="-log10(pvalue)")+
  theme(legend.title=element_text(size=14), legend.text = element_text(size=14))+
  theme(axis.title.y = element_text(margin = margin(r = 50)),axis.title.x = element_text(margin = margin(t = 20)))+
  theme(axis.text.x = element_text(face ="bold",color="black",angle=0,vjust=1))+mytheme

ggsave("SCT_Neu_pseudotime_GO_display_group.pdf",width = 11,height = 8)

BEAM_res <- BEAM(cds, branch_point = 1, cores = 1,progenitor_method = "duplicate")
BEAM_res <- BEAM_res[order(BEAM_res$qval),]
BEAM_res <- BEAM_res[,c("gene_short_name", "pval", "qval")]
BEAM_res_gene <- subset(BEAM_res,qval < 0.01)

write.csv(BEAM_res_gene,"SCT_Neu_pseudotime_BEAM_sig_qvalue_0.01gene.csv")


plot_genes_branched_heatmap(cds[BEAM_res_gene$gene_short_name,],
                            branch_point = 1,
                            num_clusters = 3,
                            cores = 1,
                            use_gene_short_name = T,
                            show_rownames = F)
BEAM_res_gene_Neu <- intersect(BEAM_res_gene$gene,Neu_marker)

plot_branch_heatmap_fixed <- function (cds_subset, 
                                       branch_point = 1, 
                                       branch_states = NULL, 
                                       branch_labels = c("Cell fate 1", "Cell fate 2"), 
                                       cluster_rows = TRUE, 
                                       hclust_method = "ward.D2", 
                                       num_clusters = 6, 
                                       hmcols = NULL, 
                                       branch_colors = c("#979797", "#F05662", "#7990C8"), 
                                       add_annotation_row = NULL, 
                                       add_annotation_col = NULL, 
                                       show_rownames = FALSE, 
                                       use_gene_short_name = TRUE, 
                                       scale_max = 3, 
                                       scale_min = -3, 
                                       norm_method = c("log", "vstExprs"), 
                                       trend_formula = "~sm.ns(Pseudotime, df=3) * Branch", 
                                       return_heatmap = FALSE, 
                                       cores = 1, 
                                       ...) 
{
  # 构建分支细胞数据集
  new_cds <- buildBranchCellDataSet(cds_subset, 
                                    branch_states = branch_states,
                                    branch_point = branch_point, 
                                    progenitor_method = "duplicate", 
                                    ...)
  new_cds@dispFitInfo <- cds_subset@dispFitInfo
  
  # 确定分支状态
  if (is.null(branch_states)) {
    progenitor_state <- subset(pData(cds_subset), Pseudotime == 0)[, "State"]
    branch_states <- setdiff(pData(cds_subset)$State, progenitor_state)
  }
  
  # 创建预测数据框架
  col_gap_ind <- 101
  newdataA <- data.frame(Pseudotime = seq(0, 100, length.out = 100), 
                         Branch = as.factor(unique(as.character(pData(new_cds)$Branch))[1]))
  newdataB <- data.frame(Pseudotime = seq(0, 100, length.out = 100), 
                         Branch = as.factor(unique(as.character(pData(new_cds)$Branch))[2]))
  
  # 生成平滑表达曲线
  BranchAB_exprs <- genSmoothCurves(new_cds[, ], 
                                    cores = cores, 
                                    trend_formula = trend_formula, 
                                    relative_expr = T, 
                                    new_data = rbind(newdataA, newdataB))
  
  # 分离两个分支的表达数据
  BranchA_exprs <- BranchAB_exprs[, 1:100]
  BranchB_exprs <- BranchAB_exprs[, 101:200]
  
  # 确定共同祖先细胞
  common_ancestor_cells <- row.names(pData(new_cds)[pData(new_cds)$State == 
                                                      setdiff(pData(new_cds)$State, branch_states), ])
  
  # 计算各分支点数
  BranchP_num <- (100 - floor(max(pData(new_cds)[common_ancestor_cells, "Pseudotime"])))
  BranchA_num <- floor(max(pData(new_cds)[common_ancestor_cells, "Pseudotime"]))
  BranchB_num <- BranchA_num
  
  # 数据标准化
  norm_method <- match.arg(norm_method)
  if (norm_method == "vstExprs") {
    BranchA_exprs <- vstExprs(new_cds, expr_matrix = BranchA_exprs)
    BranchB_exprs <- vstExprs(new_cds, expr_matrix = BranchB_exprs)
  } else if (norm_method == "log") {
    BranchA_exprs <- log10(BranchA_exprs + 1)
    BranchB_exprs <- log10(BranchB_exprs + 1)
  }
  
  # 构建热图矩阵
  heatmap_matrix <- cbind(BranchA_exprs[, (col_gap_ind - 1):1], BranchB_exprs)
  heatmap_matrix = heatmap_matrix[!apply(heatmap_matrix, 1, sd) == 0, ]
  heatmap_matrix = Matrix::t(scale(Matrix::t(heatmap_matrix), center = TRUE))
  heatmap_matrix = heatmap_matrix[is.na(row.names(heatmap_matrix)) == FALSE, ]
  heatmap_matrix[is.nan(heatmap_matrix)] = 0
  heatmap_matrix[heatmap_matrix > scale_max] = scale_max
  heatmap_matrix[heatmap_matrix < scale_min] = scale_min
  
  # 保存原始矩阵并过滤
  heatmap_matrix_ori <- heatmap_matrix
  heatmap_matrix <- heatmap_matrix[is.finite(heatmap_matrix[, 1]) & 
                                     is.finite(heatmap_matrix[, col_gap_ind]), ]
  
  # 计算行距离
  row_dist <- as.dist((1 - cor(Matrix::t(heatmap_matrix)))/2)
  row_dist[is.na(row_dist)] <- 1
  
  # 设置颜色断点
  exp_rng <- range(heatmap_matrix)
  bks <- seq(exp_rng[1] - 0.1, exp_rng[2] + 0.1, by = 0.1)
  
  # 设置默认颜色
  if (is.null(hmcols)) {
    hmcols <- colorRampPalette(c("navy", "white", "firebrick3"))(length(bks) - 1)
  }
  
  # 手动进行聚类
  if (cluster_rows) {
    row_cluster <- hclust(row_dist, method = hclust_method)
    clusters <- cutree(row_cluster, num_clusters)
  } else {
    clusters <- rep(1, nrow(heatmap_matrix))
  }
  
  # 创建行注释
  annotation_row <- data.frame(Cluster = factor(clusters))
  row.names(annotation_row) <- row.names(heatmap_matrix)
  
  if (!is.null(add_annotation_row)) {
    annotation_row <- cbind(annotation_row, 
                            add_annotation_row[row.names(annotation_row), , drop = FALSE])
  }
  
  # 创建列注释
  colnames(heatmap_matrix) <- c(1:ncol(heatmap_matrix))
  annotation_col <- data.frame(row.names = c(1:ncol(heatmap_matrix)), 
                               `Cell Type` = c(rep(branch_labels[1], BranchA_num),
                                               rep("Pre-branch", 2 * BranchP_num),
                                               rep(branch_labels[2], BranchB_num)))
  colnames(annotation_col) <- "Cell Type"
  
  if (!is.null(add_annotation_col)) {
    annotation_col <- cbind(annotation_col, 
                            add_annotation_col[row.names(annotation_col), , drop = FALSE])
  }
  
  # 设置注释颜色
  names(branch_colors) <- c("Pre-branch", branch_labels[1], branch_labels[2])
  annotation_colors = list(`Cell Type` = branch_colors)
  
  # 设置行标签
  if (use_gene_short_name == TRUE) {
    if (is.null(fData(cds_subset)$gene_short_name) == FALSE) {
      feature_label <- as.character(fData(cds_subset)[row.names(heatmap_matrix), "gene_short_name"])
      feature_label[is.na(feature_label)] <- row.names(heatmap_matrix)
    } else {
      feature_label <- row.names(heatmap_matrix)
    }
  } else {
    feature_label <- row.names(heatmap_matrix)
  }
  
  row.names(heatmap_matrix) <- feature_label
  
  # 绘制热图
  ph_res <- pheatmap(heatmap_matrix, 
                     cluster_cols = FALSE, 
                     cluster_rows = cluster_rows, 
                     show_rownames = show_rownames, 
                     show_colnames = FALSE, 
                     clustering_distance_rows = row_dist, 
                     clustering_method = hclust_method, 
                     cutree_rows = num_clusters, 
                     annotation_row = annotation_row, 
                     annotation_col = annotation_col, 
                     annotation_colors = annotation_colors, 
                     gaps_col = col_gap_ind, 
                     treeheight_row = 20, 
                     breaks = bks, 
                     fontsize = 6, 
                     color = hmcols, 
                     border_color = NA)
  
  # 返回结果
  if (return_heatmap) {
    return(list(
      heatmap = ph_res,
      clusters = clusters,  # 这是之前定义的聚类结果
      heatmap_matrix = heatmap_matrix,
      annotation_row = annotation_row
    ))
  }
}


# 使用修复后的函数
p <- plot_branch_heatmap_fixed(cds[BEAM_res_gene$gene_short_name,],
                               branch_point = 1,
                               num_clusters = 3,
                               cores = 1,
                               show_rownames = TRUE,
                               return_heatmap = T,
                               use_gene_short_name = TRUE)


clustering <- data.frame(p$annotation_row)
clustering[,1] <- as.character(clustering[,1])
colnames(clustering) <- "Gene_Clusters"
clustering$gene <- rownames(clustering)

write.csv(clustering,"SCT_pseu_BEAM_qvalue_0.01_cluster.csv")
rownames(clustering) <- clustering$gene

BEAM_res_gene_clustering <- merge(BEAM_res_gene,clustering,by="gene")

BEAM_res_gene_clustering_new <- BEAM_res_gene_clustering[,c(1,3,4,6)]
BEAM_res_gene_clustering_new2 <- BEAM_res_gene_clustering_new %>%
  left_join(cds_DGT_pseudotimegenes_sig %>% 
     select(gene_short_name, num_cells_expressed),
    by = c("gene" = "gene_short_name"))
write.csv(BEAM_res_gene_clustering_new2,"SCT_pseu_BEAM_qvalue_0.01_cluster.csv")

##单独的基因可视化

for (gene in cluster_Adult_most$gene){
  p=plot_cell_trajectory(cds, markers = gene ,use_color_gradient = T,cell_size = 1)+
    tidydr::theme_dr(xlength = 0.2, ylength = 0.2,arrow = arrow(length = unit(0.2, "inches"),type = "closed")) +
    theme(panel.grid = element_blank(),
          axis.title = element_text(face = 2,hjust = 0.03))
  ggsave(paste0(gene,"_Neu_pseu.pdf",sep=""),p,width = 6,height = 4)
}

for (gene in cluster_E_most$gene){
  p=plot_cell_trajectory(cds, markers = gene ,use_color_gradient = T,cell_size = 1)+
    tidydr::theme_dr(xlength = 0.2, ylength = 0.2,arrow = arrow(length = unit(0.2, "inches"),type = "closed")) +
    theme(panel.grid = element_blank(),
          axis.title = element_text(face = 2,hjust = 0.03))
  ggsave(paste0(gene,"_Neu_pseu.pdf",sep=""),p,width = 6,height = 4)
}


features = c("Snap25","Spp1","Calb2","Nefh","Cend1","Meg3","Fgf10","Fgf12","Ina",
             "Fabp7","Tubb2b","Prph")
pdf("SCT_Neu_pseu_cluster_gene_pseudotime.pdf",width = 8,height = 5.5)
plot_genes_in_pseudotime(cds[features,],
                         color_by = "group",
                         panel_order =c("Snap25","Spp1","Calb2","Nefh","Cend1","Meg3","Fgf10","Fgf12",
                                        "Fabp7","Tubb2b","Prph","Ina"),
                         ncol = 4,
                         cell_size = 1,
                         nrow = NULL)+scale_color_manual(values = c("#79B494", "#D67E56", "#848CBD"))
dev.off()
pdf("SCT_Neu_pseu_cluster_gene_branched_pseudotime.pdf",width = 8,height = 5.5)

plot_genes_branched_pseudotime(cds[features,],
                               branch_point = 1,
                               color_by = "group",
                               panel_order =c("Snap25","Spp1","Calb2","Nefh",
                                              "Cend1","Meg3","Fgf10","Fgf12",
                                              "Fabp7","Tubb2b","Prph","Ina"),
                               method = 'loess',
                               cell_size=2,
                               ncol = 4)+scale_color_manual(values = c("navy", "firebrick3","#79B494", "#D67E56", "#848CBD"))

dev.off()



DefaultAssay(Neu) <- "Spatial"
Idents(Neu)<-Neu$group
Idents(Neu) <- factor(Idents(Neu), levels = c("E17.5","P8","Adult"))
DotPlot(Neu, features = c("Fabp7","Tubb2b","Tubb3","Prph","Dcx","Ina","Spock3",
                          "Snap25","Nefl","Spp1","Calb2","Calb1",
                          "Nefm","Nefh","Meg3","Fgf10"))+
  RotatedAxis()+coord_flip()
ggsave("SCT_Neu_pseu_gene_dotplot.pdf",width = 5,height = 5)

rownames(BEAM_res_gene_clustering_new2) <- BEAM_res_gene_clustering_new2$gene
BEAM_res_gene_noMT <- BEAM_res_gene_clustering_new2[!grepl("^mt-", rownames(BEAM_res_gene_clustering_new2)), ]

BEAM_res_gene_noMT_main <- BEAM_res_gene_noMT[!grepl("^Gm", rownames(BEAM_res_gene_noMT)), ]
BEAM_res_gene_noMT_main <- BEAM_res_gene_noMT_main[!grepl("^Rps", rownames(BEAM_res_gene_noMT_main)), ]
BEAM_res_gene_noMT_main <- BEAM_res_gene_noMT_main[!grepl("^Rpl", rownames(BEAM_res_gene_noMT_main)), ]
BEAM_res_gene_noMT_main <- BEAM_res_gene_noMT_main[!grepl("Rik$", rownames(BEAM_res_gene_noMT_main)), ]

BEAM_res_gene_noMT_main_2 <- subset(BEAM_res_gene_noMT_main,num_cells_expressed>500)

BEAM_res_gene_noMT_main_2_cluster1 <- BEAM_res_gene_noMT_main_2[BEAM_res_gene_noMT_main_2$Gene_Clusters=="1",]
p <- plot_branch_heatmap_fixed(cds[BEAM_res_gene_noMT_main_2_cluster1$gene,],
                               branch_point = 1,
                               num_clusters = 20,
                               cores = 1,
                               show_rownames = TRUE,
                               return_heatmap = T,
                               use_gene_short_name = TRUE)


clustering_1 <- data.frame(p$annotation_row)
clustering_1[,1] <- as.character(clustering_1[,1])
colnames(clustering_1) <- "Gene_Clusters"
clustering_1$gene <- rownames(clustering_1)
clustering_1_remain <- subset(clustering_1, Gene_Clusters %in% c("6", "2", "10","13"))


BEAM_res_gene_noMT_main_2_cluster3 <- BEAM_res_gene_noMT_main_2[BEAM_res_gene_noMT_main_2$Gene_Clusters=="3",]
p <- plot_branch_heatmap_fixed(cds[BEAM_res_gene_noMT_main_2_cluster3$gene,],
                               branch_point = 1,
                               num_clusters = 2,
                               cores = 1,
                               show_rownames = TRUE,
                               return_heatmap = T,
                               use_gene_short_name = TRUE)


clustering_3 <- data.frame(p$annotation_row)
clustering_3[,1] <- as.character(clustering_3[,1])
colnames(clustering_3) <- "Gene_Clusters"
clustering_3$gene <- rownames(clustering_3)
clustering_3_remain <- clustering_3[clustering_3$Gene_Clusters=="1",]

BEAM_res_gene_noMT_main_2_cluster2 <- BEAM_res_gene_noMT_main_2[BEAM_res_gene_noMT_main_2$Gene_Clusters=="2",]

BEAM_res_gene_noMT_main_2_select <- union(clustering_1_remain$gene,BEAM_res_gene_noMT_main_2_cluster2$gene)
BEAM_res_gene_noMT_main_2_select <- union(BEAM_res_gene_noMT_main_2_select,clustering_3_remain$gene)

BEAM_res_gene_noMT_main_2_select_info <- BEAM_res_gene_noMT[BEAM_res_gene_noMT$gene%in%BEAM_res_gene_noMT_main_2_select,]
write.csv(BEAM_res_gene_noMT_main_2_select_info,"SCT_Neu_BEAM_heatmap_gene_info.csv")

plot_genes_branched_heatmap(cds[BEAM_res_gene_noMT_main_2_select,],
                            branch_point = 1,
                            num_clusters = 3,
                            cores = 1,
                            use_gene_short_name = T,
                            show_rownames = T)
p <- plot_branch_heatmap_fixed(cds[BEAM_res_gene_noMT_main_2_select,],
                               branch_point = 1,
                               num_clusters = 3,
                               cores = 1,
                               show_rownames = TRUE,
                               return_heatmap = T,
                               use_gene_short_name = TRUE)


clustering <- data.frame(p$annotation_row)
clustering[,1] <- as.character(clustering[,1])
colnames(clustering) <- "Gene_Clusters"
clustering$gene <- rownames(clustering)

library(clusterProfiler)
library(msigdbr)
library(GSVA)
library(org.Mm.eg.db)
library(enrichR)

GO_all<-data.frame(Description=character(),
                   qvalue=numeric(),
                   pvalue=numeric(),
                   RichFactor=numeric(),
                   Count=numeric(),
                   Ontology=character(),
                   gene=character(),
                   group=character())
for (i in unique(clustering$Gene_Clusters)) {
  gene=clustering[clustering$Gene_Clusters==i,]$gene
  gene_entrez_id2GO <- clusterProfiler::bitr(gene, "SYMBOL", "ENTREZID", "org.Mm.eg.db", drop = TRUE)$ENTREZID    
  eG <- enrichGO(gene = gene_entrez_id2GO, #需要分析的基因的EntrezID
                 OrgDb = org.Mm.eg.db, #人基因数据库
                 pvalueCutoff =0.01, #设置pvalue界值
                 qvalueCutoff = 0.01, #设置qvalue界值(FDR校正后的p值）
                 ont="all", #选择功能富集的类型，可选BP、MF、CC，这里选择all。
                 readable =T)
  GO <- data.frame(eG@result)
  GO$group <- rep(i,nrow(GO))
  GO_new <- data.frame(Description=GO$Description,
                       qvalue=GO$qvalue,
                       pvalue=GO$pvalue,
                       RichFactor=GO$RichFactor,
                       Count=GO$Count,
                       Ontology=GO$ONTOLOGY,
                       gene=GO$geneID,
                       group=GO$group)
  GO_all <- rbind(GO_all,GO_new)
}
write.csv(GO_all,"SCT_Neu_BEAM_GO_results.csv")
GO_terms <- GO_all%>%group_by(group)%>%top_n(10,-pvalue)

ggplot(GO_terms)+
  geom_bar(aes(-log10(as.numeric(pvalue)),fct_reorder(Description,group),fill = Ontology),stat="identity") + #将Term按照Ontology排序
  theme_bw()+
  theme(text=element_text(family='serif',size =10,face = 'bold'),axis.text = element_text(family='serif',size =10,face = 'bold',colour = as.factor(group))) + 
  #geom_text(data=subset(GO_terms,qvalue<0.05,c('Count','Description')),aes(x=Count+3,y=as.factor(Description),label='*')) + #对于qvalue<0.05的term添加*
  labs(x='-log10(pvalue)',y='GO term') + 
  scale_fill_manual(values=c(BP = "#79B494", CC = "#D67E56", MF = "#848CBD"))

ggplot(GO_terms) +
  geom_bar(aes(-log10(as.numeric(pvalue)), fct_reorder(Description, group), fill = Ontology), stat = "identity") +
  theme_bw() +
  theme(
    text = element_text(family = 'serif', size = 10, face = 'bold'),
    axis.text.y = element_text(family = 'serif', size = 10, face = 'bold', 
                               colour = as.character(GO_terms$group)),  # 按分组设置y轴标签颜色
    axis.text.x = element_text(family = 'serif', size = 10, face = 'bold', colour = "black")
  ) +
  labs(x = '-log10(pvalue)', y = 'GO term') +
  scale_fill_manual(values = c(BP = "#79B494", CC = "#D67E56", MF = "#848CBD"))
ggsave("SCT_neu_BEAM_top10_GO_barplot.pdf",width = 8,height = 8)














####20250920
Neu <- readRDS("SCT_group_neu.Rds")
Idents(Neu)<-Neu$group

integrated_data <- SeuratObject::LayerData(Neu, assay = "integrated")
integrated_cells <- colnames(integrated_data)
p_data <- Neu@meta.data[integrated_cells, ]
f_data <- data.frame(gene_short_name = rownames(integrated_data),
                     row.names = rownames(integrated_data))
pd<-new('AnnotatedDataFrame',data=p_data)
fd<-new('AnnotatedDataFrame',data=f_data)
##创建对象
cds<-newCellDataSet(integrated_data,
                    phenoData = pd,
                    featureData = fd,
                    lowerDetectionLimit = 0.5,
                    expressionFamily = negbinomial.size())
##归一化处理
cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)
##数据过滤
cds<-detectGenes(cds,min_expr = 0.1)
##选择过程基因
##选择seurat高变异的基因
expression_gene <- VariableFeatures(Neu)
cds <- setOrderingFilter(cds,expression_gene)
plot_ordering_genes(cds)
##选择COSG所计算的各分组高表达的基因

##
cds <- setOrderingFilter(cds,Neu_marker)
plot_ordering_genes(cds)
##

##选择monocle计算的差异基因
diff_test_res <- differentialGeneTest(cds, fullModelFormulaStr = "~group",cores = 4)
ordering_genes <- rownames(subset(diff_test_res, qval < 0.05))

cds <- setOrderingFilter(cds,ordering_genes)
plot_ordering_genes(cds)

##降维

cds <- reduceDimension(cds,method='DDRTree')
##排序
cds <- orderCells(cds)
# 可视化轨迹，按伪时间着色
plot_cell_trajectory(cds, color_by = "Pseudotime")

# 可视化轨迹，按样本着色
plot_cell_trajectory(cds, color_by = "group")

# 可视化轨迹，按样本着色
plot_cell_trajectory(cds, color_by = "State")

##热图绘制
cds_DGT_pseudotimegenes <- differentialGeneTest(cds, fullModelFormulaStr = "~group")
cds_DGT_pseudotimegenes_sig <- subset(cds_DGT_pseudotimegenes, qval < 0.05)
cds_DGT_pseudotimegenes_sig <- subset(cds_DGT_pseudotimegenes, pval < 0.05)
cds_DGT_pseudotimegenes_sig_GO_gene <- cds_DGT_pseudotimegenes[remain_gene_new,]

write.csv(cds_DGT_pseudotimegenes_sig,"SCT_neu_peudotimegenes_qvalue_0.05.csv")
write.csv(cds_DGT_pseudotimegenes_sig,"SCT_neu_peudotimegenes_pvalue_0.05.csv")
write.csv(cds_DGT_pseudotimegenes_sig_GO_gene,"SCT_neu_peudotimegenes_GO_gene.csv")

cds_DGT_pseudotimegenes_sig <- subset(cds_DGT_pseudotimegenes, pval < 0.05)
remain_gene_new

# 然后绘制热图
plot_pseudotime_heatmap(cds[cds_DGT_pseudotimegenes_sig$gene_short_name,], 
                        num_cluster = 3, 
                        show_rownames = F, 
                        return_heatmap = T,
                        #add_annotation_col = annotation_col,
                        hmcols = colorRampPalette(c("navy","white","firebrick3"))(100))

marker <- FindAllMarkers(Neu, only.pos = T,logfc.threshold = 0.5)
marker <- marker[which(marker$p_val_adj<0.05),]
top15 <- marker %>% group_by(cluster) %>% top_n(n = 100, wt = avg_log2FC)
top15_ordergene <- cds_DGT_pseudotimegenes_sig[top15$gene, ]
top15_ordergene <- top15_ordergene[!grepl("^NA", rownames(top15_ordergene)), ]
Time_genes <- top15_ordergene %>% pull(gene_short_name) %>% as.character()
Time_genes <- unique(Time_genes)

plot_pseudotime_heatmap(cds[Time_genes,], 
                        num_cluster = 3, 
                        show_rownames = T, 
                        return_heatmap = T,
                        hmcols = colorRampPalette(c("navy","white","firebrick3"))(100))


top15_ordergene <- cds_DGT_pseudotimegenes_sig[top15$gene, ]
top15_ordergene <- top15_ordergene[!grepl("^NA", rownames(top15_ordergene)), ]
Time_genes <- top15_ordergene %>% pull(gene_short_name) %>% as.character()
Time_genes <- unique(Time_genes)
add_gene <- c("Snap25","Tuba1a","Tubb2b","Spp1","Calb2")
select_gene <- union(add_gene,Time_genes)


plot_pseudotime_heatmap(cds[select_gene,], 
                        num_cluster = 3, 
                        show_rownames = T, 
                        return_heatmap = T,
                        hmcols = colorRampPalette(c("navy","white","firebrick3"))(100))

COSG_ordergene <- cds_DGT_pseudotimegenes_sig[Neu_marker, ]
COSG_ordergene <- top15_ordergene[!grepl("^NA", rownames(top15_ordergene)), ]
Time_genes <- top15_ordergene %>% pull(gene_short_name) %>% as.character()
Time_genes <- unique(Time_genes)


plot_pseudotime_heatmap(cds[Time_genes,], 
                        num_cluster = 3, 
                        show_rownames = T, 
                        return_heatmap = F,
                        hmcols = colorRampPalette(c("navy","white","firebrick3"))(100))

time_genes <- cds_DGT_pseudotimegenes_sig$gene_short_name


##GO富集分析
library(clusterProfiler)
library(msigdbr)
library(GSVA)
library(org.Mm.eg.db)
library(enrichR)
View(msigdbr_collections())
setwd("../")
## 获取各分支基因
cds_DGT_pseudotimegenes_sig <- subset(cds_DGT_pseudotimegenes, pval < 0.05)

gene=cds_DGT_pseudotimegenes_sig$gene_short_name

## 把SYMBOL改为ENTREZID
gene_entrez_id2GO <- clusterProfiler::bitr(gene, "SYMBOL", "ENTREZID", "org.Mm.eg.db", drop = TRUE)$ENTREZID    

## GO
eG <- enrichGO(gene = gene_entrez_id2GO, #需要分析的基因的EntrezID
               OrgDb = org.Mm.eg.db, #人基因数据库
               pvalueCutoff =0.01, #设置pvalue界值
               qvalueCutoff = 0.01, #设置qvalue界值(FDR校正后的p值）
               ont="all", #选择功能富集的类型，可选BP、MF、CC，这里选择all。
               readable =T)
barplot(eG, x = "GeneRatio", color = "p.adjust", #默认参数（x和color可以根据eG里面的内容更改）
          showCategory =10, #只显示前10
          split="ONTOLOGY") + #以ONTOLOGY类型分开
    facet_grid(ONTOLOGY~., scale='free') +labs(title = "eGO_barplot_Neu_cluster1_top10GO")
Neu_pseu_GO <- data.frame(eG@result)

#结果优化
library(ggplot2)
library(forcats)
GO <- data.frame(Term=Neu_pseu_GO$Description,
                 qvalue=Neu_pseu_GO$qvalue,
                 Count=Neu_pseu_GO$Count,
                 Ontology=Neu_pseu_GO$ONTOLOGY,
                 gene=Neu_pseu_GO$geneID)
write.csv(GO,"SCT_Neu_pseu_sig_gene_GO_select.csv")
GO <- read.csv("SCT_Neu_pseu_sig_gene_GO_select.csv")

ggplot(GO)+
  geom_bar(aes(Count,fct_reorder(Term,Ontology),fill = Ontology),stat="identity") + #将Term按照Ontology排序
  theme_bw()+
  theme(text=element_text(family='serif',size =10,face = 'bold'),axis.text = element_text(family='serif',size =10,face = 'bold',colour = 'black')) + 
  geom_text(data=subset(GO,qvalue<0.05,c('Count','Term')),aes(x=Count+3,y=as.factor(Term),label='*')) + #对于qvalue<0.05的term添加*
  labs(x='Number of genes',y='GO term') + 
  scale_fill_manual(values=c(BP = "#79B494", CC = "#D67E56", MF = "#848CBD"))

DefaultAssay(Neu)<-"integrated"
DotPlot(Neu, features = c("Snap25","Nefl","Spp1","Calb2","Rps7",
                          "Calb1","Fgf14","Gjb2","Pllp","Gjb6","Epyc",
                          "Cplx1","Syt7","Syt2","Vamp1","S100b","Ptgds"))+ RotatedAxis()+
  theme(legend.position = "none")+
  labs(title = "SCT_Neu_pseu_expression")
ggsave("SCT_Neu_pseu_expression.pdf")
Neu_pseu_GO <- data.frame(eG@result)
write.csv(Neu_pseu_GO,"SCT_Neu_pseu_sig_gene_GO.csv")





Neu_pseu_GO_gene <- unique(unlist(strsplit(Neu_pseu_GO$geneID, "/")))

pdf("SCT_Neu_pseudotime_heatmap.pdf",height = 6,width = 4)
##筛选可视化基因
plot_pseudotime_heatmap(cds[Neu_pseu_GO_gene,], 
                        num_cluster = 10, 
                        show_rownames = F, 
                        return_heatmap = T,
                        hmcols = colorRampPalette(c("navy","white","firebrick3"))(100))
dev.off()
b=plot_pseudotime_heatmap(cds[Neu_pseu_GO_gene,], 
                          num_cluster = 10, 
                          show_rownames = F, 
                          return_heatmap = T,
                          hmcols = colorRampPalette(c("navy","white","firebrick3"))(100))
clusters_all <- cutree(b$tree_row, k = 10)
clustering_all <- data.frame(clusters_all)
clustering_all[,1] <- as.character(clustering_all[,1])
colnames(clustering_all) <- "Gene_Clusters"
clustering_all$gene <- rownames(clustering_all)


Neu_cluster_gene <- clustering_all[clustering_all$Gene_Clusters=="10",]$gene
plot_pseudotime_heatmap(cds[remain_gene_new,], 
                          num_cluster = 3, 
                          show_rownames = T, 
                          return_heatmap = T,
                          hmcols = colorRampPalette(c("navy","white","firebrick3"))(100))
p=plot_pseudotime_heatmap(cds[remain_gene,], 
                          num_cluster = 6, 
                          show_rownames = T, 
                          return_heatmap = T,
                          hmcols = colorRampPalette(c("navy","white","firebrick3"))(100))

clusters <- cutree(p$tree_row, k = 6)
clustering <- data.frame(clusters)
clustering[,1] <- as.character(clustering[,1])
colnames(clustering) <- "Gene_Clusters"
clustering$gene <- rownames(clustering)

remain_gene <- subset(clustering,Gene_Clusters!="1")$gene
remain_gene_new <- subset(clustering,Gene_Clusters!="5")$gene
remain_gene_new <- subset(clustering_all,Gene_Clusters=="10")$gene

remain_gene <- union(remain_gene_new,remain_gene)
plot_pseudotime_heatmap(cds[remain_gene,], 
                        num_cluster = 3, 
                        show_rownames = T, 
                        return_heatmap = T,
                        hmcols = colorRampPalette(c("navy","white","firebrick3"))(100))


##GO富集分析
library(clusterProfiler)
library(msigdbr)
library(GSVA)
library(org.Mm.eg.db)
library(enrichR)
View(msigdbr_collections())
setwd("../")
## 获取各分支基因
p=plot_pseudotime_heatmap(cds[remain_gene_new,], 
                          num_cluster = 3, 
                          show_rownames = T, 
                          return_heatmap = T,
                          hmcols = colorRampPalette(c("navy","white","firebrick3"))(100))

clusters <- cutree(p$tree_row, k = 3)
clustering <- data.frame(clusters)
clustering[,1] <- as.character(clustering[,1])
colnames(clustering) <- "Gene_Clusters"
clustering$gene <- rownames(clustering)

GO_all<-data.frame(Description=character(),
                   qvalue=numeric(),
                   pvalue=numeric(),
                   Count=numeric(),
                   ONTOLOGY=character(),
                   geneID=character(),
                   group=character())
for (i in 1:3) {
  gene=clustering[clustering$Gene_Clusters==i,]$gene
  gene_entrez_id2GO <- clusterProfiler::bitr(gene, "SYMBOL", "ENTREZID", "org.Mm.eg.db", drop = TRUE)$ENTREZID    
  eG <- enrichGO(gene = gene_entrez_id2GO, #需要分析的基因的EntrezID
                 OrgDb = org.Mm.eg.db, #人基因数据库
                 pvalueCutoff =0.01, #设置pvalue界值
                 qvalueCutoff = 0.01, #设置qvalue界值(FDR校正后的p值）
                 ont="all", #选择功能富集的类型，可选BP、MF、CC，这里选择all。
                 readable =T)
  GO <- data.frame(eG@result)
  GO$group <- rep(i,nrow(GO))
  GO_new <- data.frame(Term=GO$Description,
                   qvalue=GO$qvalue,
                   pvalue=GO$pvalue,
                   Count=GO$Count,
                   Ontology=GO$ONTOLOGY,
                   gene=GO$geneID,
                   cluster=GO$group)
  GO_all <- rbind(GO_all,GO_new)
}
write.csv(GO_all,"SCT_Neu_pseu_sig_gene_GO_select_cluster.csv")
GO_all <- read.csv("SCT_Neu_pseu_sig_gene_GO_select_cluster.csv")

library(ggplot2)
library(forcats)

test <- GO[GO$cluster==1,]
p1<-ggplot(test)+
  geom_bar(aes(Count,fct_reorder(Term,Ontology),fill = Ontology),stat="identity") + #将Term按照Ontology排序
  theme_bw()+
  theme(text=element_text(family='serif',size =10,face = 'bold'),
        axis.text = element_text(family='serif',size =10,face = 'bold',colour = 'black'),
        panel.grid = element_blank()) + 
  #geom_text(data=subset(test,qvalue<0.05,c('Count','Term')),aes(x=Count+3,y=as.factor(Term),label='*')) + #对于qvalue<0.05的term添加*
  labs(x='Number of genes',y='GO term',title = "SCT_Neu_cluster1_gene_GO_term") + 
  scale_fill_manual(values=c(BP = "#79B494", CC = "#D67E56", MF = "#848CBD"))
test <- GO[GO$cluster==2,]
p2<-ggplot(test)+
  geom_bar(aes(Count,fct_reorder(Term,Ontology),fill = Ontology),stat="identity") + #将Term按照Ontology排序
  theme_bw()+
  theme(text=element_text(family='serif',size =10,face = 'bold'),
        axis.text = element_text(family='serif',size =10,face = 'bold',colour = 'black'),
        panel.grid = element_blank()) + 
  #geom_text(data=subset(test,qvalue<0.05,c('Count','Term')),aes(x=Count+3,y=as.factor(Term),label='*')) + #对于qvalue<0.05的term添加*
  labs(x='Number of genes',y='GO term',title = "SCT_Neu_cluster2_gene_GO_term") + 
  scale_fill_manual(values=c(BP = "#79B494", CC = "#D67E56", MF = "#848CBD"))
test <- GO[GO$cluster==3,]
p3<-ggplot(test)+
  geom_bar(aes(Count,fct_reorder(Term,Ontology),fill = Ontology),stat="identity") + #将Term按照Ontology排序
  theme_bw()+
  theme(text=element_text(family='serif',size =10,face = 'bold'),
        axis.text = element_text(family='serif',size =10,face = 'bold',colour = 'black'),
        panel.grid = element_blank()) + 
  #geom_text(data=subset(test,qvalue<0.05,c('Count','Term')),aes(x=Count+3,y=as.factor(Term),label='*')) + #对于qvalue<0.05的term添加*
  labs(x='Number of genes',y='GO term',title = "SCT_Neu_cluster3_gene_GO_term") + 
  scale_fill_manual(values=c(BP = "#79B494", CC = "#D67E56", MF = "#848CBD"))
cowplot::plot_grid(p1,p2,p3,ncol = 1,align = 'hv')
ggsave("SCT_Neu_pseu_cluster_GO_term_log10pvalue.pdf",width = 8,height = 8)


##热图数据提取
b=plot_pseudotime_heatmap(cds[Neu_pseu_GO_gene,], 
                        num_cluster = 3, 
                        show_rownames = T, 
                        return_heatmap = T,
                        hmcols = colorRampPalette(c("navy","white","firebrick3"))(100))

newdata<-data.frame(Pseudotime = seq(min(pData(cds)$Pseudotime), 
                                     max(pData(cds)$Pseudotime), length.out = 100))
genSmoothCurves_mat<-genSmoothCurves(cds[Neu_pseu_GO_gene,],
                                     new_data = newdata,
                                     cores = 10)
pheatmap(log10(genSmoothCurves_mat[b$tree_row$order,]+1),
         scale = "row",
         cluster_rows = F,
         breaks = seq(-3,3,6/100),
         cluster_cols = F,
         #annotation_col = annocol,
         #annotation_colors = annocolor,
         show_rownames = F,
         show_colnames = F,
         color = colorRampPalette(c("navy","white","firebrick3"))(100))



#Monocle2提供了一种特殊的统计检验方法：branched expression analysis modeling（BEAM），可以对不同的分支事件进行分析。
BEAM_res <- BEAM(cds, branch_point = 1, cores = 1,progenitor_method = "duplicate")
BEAM_res <- BEAM_res[order(BEAM_res$qval),]
BEAM_res <- BEAM_res[,c("gene_short_name", "pval", "qval")]
BEAM_res_gene <- subset(BEAM_res,qval < 0.05)
pdf("heatmap_1.pdf",height = 6,width = 4)
plot_genes_branched_heatmap(cds[Neu_pseu_GO_gene,],
                            branch_point = 1,
                            num_clusters = 3,
                            cores = 1,
                            use_gene_short_name = T,
                            show_rownames = T,
                            hmcols = colorRampPalette(c("navy","white","firebrick3"))(100))
dev.off()
# 保存整个cds对象（包含所有分析结果）
saveRDS(cds, file = "SCT_Neu_pseudotime_analysis_results_2.Rds")

##结果可视化部分
cds <- readRDS("SCT_Neu_pseudotime_analysis_results.Rds")

# 可视化轨迹，按伪时间着色
plot_cell_trajectory(cds, color_by = "Pseudotime")+
  tidydr::theme_dr(xlength = 0.2, ylength = 0.2,arrow = arrow(length = unit(0.2, "inches"),type = "closed")) +
  theme(panel.grid = element_blank(),
        axis.title = element_text(face = 2,hjust = 0.03))
ggsave("SCT_Neu_pseu_Pseudotime.pdf",width = 5.5,height = 5)

# 可视化轨迹，按样本着色
plot_cell_trajectory(cds, color_by = "group")+scale_colour_manual(values = c("#79B494", "#D67E56",  "#848CBD"))+
  tidydr::theme_dr(xlength = 0.2, ylength = 0.2,arrow = arrow(length = unit(0.2, "inches"),type = "closed")) +
  theme(panel.grid = element_blank(),
        axis.title = element_text(face = 2,hjust = 0.03))
ggsave("SCT_Neu_pseu_group.pdf",width = 5.5,height = 5)
# 可视化轨迹，按样本着色
plot_cell_trajectory(cds, color_by = "State")+
  tidydr::theme_dr(xlength = 0.2, ylength = 0.2,arrow = arrow(length = unit(0.2, "inches"),type = "closed")) +
  theme(panel.grid = element_blank(),
        axis.title = element_text(face = 2,hjust = 0.03))&NoLegend()
ggsave("SCT_Neu_pseu_State.pdf",width = 5.5,height = 5)

DefaultAssay(Neu)<-"integrated"
DotPlot(Neu, features = c("Snap25","Nefl","Spp1","Calb2","Tubb2b",
                          
                          "Cplx1","Fgf1","Syt7","Vamp1","S100b","Ptgds"))+
  RotatedAxis()+coord_flip()
ggsave("SCT_Neu_pseu_gene_dotplot.pdf",width = 5,height = 6)



##剔除较远的细胞

Neu_idents <- data.frame(Neu@reductions$umap@cell.embeddings)
cells_to_exclude <- rownames(Neu_idents[Neu_idents$umap_1>-4,])
Neu_filtered <- Neu[, !colnames(Neu) %in% cells_to_exclude]
DimPlot(Neu_filtered, reduction = "umap",label = F,group.by = "group")

Neu_filtered <- ScaleData(Neu_filtered, verbose = FALSE)
Neu_filtered <- RunPCA(Neu_filtered, npcs = 50, verbose = FALSE)
Neu_filtered <- RunUMAP(Neu_filtered, reduction = "pca", dims = 1:50)
Neu_filtered <- FindNeighbors(Neu_filtered, reduction = "pca", dims = 1:50)
Neu_filtered <- FindClusters(Neu_filtered, resolution = 0.1,verbose = FALSE)
DimPlot(Neu_filtered, reduction = "umap", group.by = c("group"))

DimPlot(Neu, reduction = "umap", group.by = c("group"))
DimPlot(Neu_filtered, reduction = "umap",pt.size = 2, cols = c("#79B494", "#D67E56",  "#848CBD"),label = F,label.size = 5,group.by = c("group"))+
  tidydr::theme_dr(xlength = 0.2, ylength = 0.2,arrow = arrow(length = unit(0.2, "inches"),type = "closed")) +
  theme(panel.grid = element_blank(),
        axis.title = element_text(face = 2,hjust = 0.03))
ggsave("SCT_Neu_group_Umap.pdf",width = 5.5,height = 5)

Neu_pseu_GO<-read.csv("SCT_Neu_pseu_sig_gene_GO.csv")
Neu_pseu_GO_gene <- unique(unlist(strsplit(Neu_pseu_GO$geneID, "/")))
pdf("SCT_Neu_pseudotime_heatmap.pdf",height = 6,width = 4)
plot_pseudotime_heatmap(cds[Neu_pseu_GO_gene,], 
                        num_cluster = 10, 
                        show_rownames = F, 
                        return_heatmap = T,
                        hmcols = colorRampPalette(c("navy","white","firebrick3"))(100))
dev.off()
dev.new()
GO<-read_xlsx("SCT_Neu_pseu_sig_gene_GO_select.xlsx")

ggplot(GO)+
  geom_bar(aes(Count,fct_reorder(Term,Ontology),fill = Ontology),stat="identity") + #将Term按照Ontology排序
  theme_bw()+
  theme(text=element_text(family='serif',size =10,face = 'bold'),axis.text = element_text(family='serif',size =10,face = 'bold',colour = 'black')) + 
  geom_text(data=subset(GO,qvalue<0.05,c('Count','Term')),aes(x=Count+3,y=as.factor(Term),label='*')) + #对于qvalue<0.05的term添加*
  labs(x='Number of genes',y='GO term') + 
  scale_fill_manual(values=c(BP = "#79B494", CC = "#D67E56", MF = "#848CBD"))





##分支节点基因热图
BEAM_res <- BEAM(cds, branch_point = 1, cores = 1,progenitor_method = "duplicate")
BEAM_res <- BEAM_res[order(BEAM_res$qval),]
BEAM_res <- BEAM_res[,c("gene_short_name", "pval", "qval")]
BEAM_res_gene <- subset(BEAM_res,pval < 0.05)
plot_genes_branched_heatmap(cds[BEAM_res_gene$gene_short_name,],
                            branch_point = 1,
                            num_clusters = 3,
                            cores = 1,
                            use_gene_short_name = T,
                            show_rownames = T)

write.csv(BEAM_res_gene,"SCT_Neu_pseudotime_BEAM_sig_gene.csv")
pdf("SCT_Neu_pseudotime_branched_heatmap.pdf",height = 6,width = 4)

plot_genes_branched_heatmap(cds[BEAM_res_gene$gene_short_name,],
                            branch_point = 1,
                            num_clusters = 3,
                            cores = 1,
                            show_rownames = T,
                            use_gene_short_name = T)
dev.off()
##分支节点基因热图GO分析
library(clusterProfiler)
library(msigdbr)
library(GSVA)
library(org.Mm.eg.db)
library(enrichR)
View(msigdbr_collections())
setwd("../")

plot_branch_heatmap_fixed <- function (cds_subset, 
                                       branch_point = 1, 
                                       branch_states = NULL, 
                                       branch_labels = c("Cell fate 1", "Cell fate 2"), 
                                       cluster_rows = TRUE, 
                                       hclust_method = "ward.D2", 
                                       num_clusters = 6, 
                                       hmcols = NULL, 
                                       branch_colors = c("#979797", "#F05662", "#7990C8"), 
                                       add_annotation_row = NULL, 
                                       add_annotation_col = NULL, 
                                       show_rownames = FALSE, 
                                       use_gene_short_name = TRUE, 
                                       scale_max = 3, 
                                       scale_min = -3, 
                                       norm_method = c("log", "vstExprs"), 
                                       trend_formula = "~sm.ns(Pseudotime, df=3) * Branch", 
                                       return_heatmap = FALSE, 
                                       cores = 1, 
                                       ...) 
{
  # 构建分支细胞数据集
  new_cds <- buildBranchCellDataSet(cds_subset, 
                                    branch_states = branch_states,
                                    branch_point = branch_point, 
                                    progenitor_method = "duplicate", 
                                    ...)
  new_cds@dispFitInfo <- cds_subset@dispFitInfo
  
  # 确定分支状态
  if (is.null(branch_states)) {
    progenitor_state <- subset(pData(cds_subset), Pseudotime == 0)[, "State"]
    branch_states <- setdiff(pData(cds_subset)$State, progenitor_state)
  }
  
  # 创建预测数据框架
  col_gap_ind <- 101
  newdataA <- data.frame(Pseudotime = seq(0, 100, length.out = 100), 
                         Branch = as.factor(unique(as.character(pData(new_cds)$Branch))[1]))
  newdataB <- data.frame(Pseudotime = seq(0, 100, length.out = 100), 
                         Branch = as.factor(unique(as.character(pData(new_cds)$Branch))[2]))
  
  # 生成平滑表达曲线
  BranchAB_exprs <- genSmoothCurves(new_cds[, ], 
                                    cores = cores, 
                                    trend_formula = trend_formula, 
                                    relative_expr = T, 
                                    new_data = rbind(newdataA, newdataB))
  
  # 分离两个分支的表达数据
  BranchA_exprs <- BranchAB_exprs[, 1:100]
  BranchB_exprs <- BranchAB_exprs[, 101:200]
  
  # 确定共同祖先细胞
  common_ancestor_cells <- row.names(pData(new_cds)[pData(new_cds)$State == 
                                                      setdiff(pData(new_cds)$State, branch_states), ])
  
  # 计算各分支点数
  BranchP_num <- (100 - floor(max(pData(new_cds)[common_ancestor_cells, "Pseudotime"])))
  BranchA_num <- floor(max(pData(new_cds)[common_ancestor_cells, "Pseudotime"]))
  BranchB_num <- BranchA_num
  
  # 数据标准化
  norm_method <- match.arg(norm_method)
  if (norm_method == "vstExprs") {
    BranchA_exprs <- vstExprs(new_cds, expr_matrix = BranchA_exprs)
    BranchB_exprs <- vstExprs(new_cds, expr_matrix = BranchB_exprs)
  } else if (norm_method == "log") {
    BranchA_exprs <- log10(BranchA_exprs + 1)
    BranchB_exprs <- log10(BranchB_exprs + 1)
  }
  
  # 构建热图矩阵
  heatmap_matrix <- cbind(BranchA_exprs[, (col_gap_ind - 1):1], BranchB_exprs)
  heatmap_matrix = heatmap_matrix[!apply(heatmap_matrix, 1, sd) == 0, ]
  heatmap_matrix = Matrix::t(scale(Matrix::t(heatmap_matrix), center = TRUE))
  heatmap_matrix = heatmap_matrix[is.na(row.names(heatmap_matrix)) == FALSE, ]
  heatmap_matrix[is.nan(heatmap_matrix)] = 0
  heatmap_matrix[heatmap_matrix > scale_max] = scale_max
  heatmap_matrix[heatmap_matrix < scale_min] = scale_min
  
  # 保存原始矩阵并过滤
  heatmap_matrix_ori <- heatmap_matrix
  heatmap_matrix <- heatmap_matrix[is.finite(heatmap_matrix[, 1]) & 
                                     is.finite(heatmap_matrix[, col_gap_ind]), ]
  
  # 计算行距离
  row_dist <- as.dist((1 - cor(Matrix::t(heatmap_matrix)))/2)
  row_dist[is.na(row_dist)] <- 1
  
  # 设置颜色断点
  exp_rng <- range(heatmap_matrix)
  bks <- seq(exp_rng[1] - 0.1, exp_rng[2] + 0.1, by = 0.1)
  
  # 设置默认颜色
  if (is.null(hmcols)) {
    hmcols <- colorRampPalette(c("navy", "white", "firebrick3"))(length(bks) - 1)
  }
  
  # 手动进行聚类
  if (cluster_rows) {
    row_cluster <- hclust(row_dist, method = hclust_method)
    clusters <- cutree(row_cluster, num_clusters)
  } else {
    clusters <- rep(1, nrow(heatmap_matrix))
  }
  
  # 创建行注释
  annotation_row <- data.frame(Cluster = factor(clusters))
  row.names(annotation_row) <- row.names(heatmap_matrix)
  
  if (!is.null(add_annotation_row)) {
    annotation_row <- cbind(annotation_row, 
                            add_annotation_row[row.names(annotation_row), , drop = FALSE])
  }
  
  # 创建列注释
  colnames(heatmap_matrix) <- c(1:ncol(heatmap_matrix))
  annotation_col <- data.frame(row.names = c(1:ncol(heatmap_matrix)), 
                               `Cell Type` = c(rep(branch_labels[1], BranchA_num),
                                               rep("Pre-branch", 2 * BranchP_num),
                                               rep(branch_labels[2], BranchB_num)))
  colnames(annotation_col) <- "Cell Type"
  
  if (!is.null(add_annotation_col)) {
    annotation_col <- cbind(annotation_col, 
                            add_annotation_col[row.names(annotation_col), , drop = FALSE])
  }
  
  # 设置注释颜色
  names(branch_colors) <- c("Pre-branch", branch_labels[1], branch_labels[2])
  annotation_colors = list(`Cell Type` = branch_colors)
  
  # 设置行标签
  if (use_gene_short_name == TRUE) {
    if (is.null(fData(cds_subset)$gene_short_name) == FALSE) {
      feature_label <- as.character(fData(cds_subset)[row.names(heatmap_matrix), "gene_short_name"])
      feature_label[is.na(feature_label)] <- row.names(heatmap_matrix)
    } else {
      feature_label <- row.names(heatmap_matrix)
    }
  } else {
    feature_label <- row.names(heatmap_matrix)
  }
  
  row.names(heatmap_matrix) <- feature_label
  
  # 绘制热图
  ph_res <- pheatmap(heatmap_matrix, 
                     cluster_cols = FALSE, 
                     cluster_rows = cluster_rows, 
                     show_rownames = show_rownames, 
                     show_colnames = FALSE, 
                     clustering_distance_rows = row_dist, 
                     clustering_method = hclust_method, 
                     cutree_rows = num_clusters, 
                     annotation_row = annotation_row, 
                     annotation_col = annotation_col, 
                     annotation_colors = annotation_colors, 
                     gaps_col = col_gap_ind, 
                     treeheight_row = 20, 
                     breaks = bks, 
                     fontsize = 6, 
                     color = hmcols, 
                     border_color = NA)
  
  # 返回结果
  if (return_heatmap) {
    return(list(
      heatmap = ph_res,
      clusters = clusters,  # 这是之前定义的聚类结果
      heatmap_matrix = heatmap_matrix,
      annotation_row = annotation_row
    ))
  }
}


# 使用修复后的函数
p <- plot_branch_heatmap_fixed(cds[BEAM_res_gene$gene_short_name,],
                               branch_point = 1,
                               num_clusters = 3,
                               cores = 1,
                               show_rownames = TRUE,
                               return_heatmap = T,
                               use_gene_short_name = TRUE)


clustering <- data.frame(p$annotation_row)
clustering[,1] <- as.character(clustering[,1])
colnames(clustering) <- "Gene_Clusters"
clustering$gene <- rownames(clustering)

GO_all<-data.frame(Description=character(),
                   qvalue=numeric(),
                   pvalue=numeric(),
                   Count=numeric(),
                   ONTOLOGY=character(),
                   geneID=character(),
                   group=character())
for (i in 1:3) {
  gene=clustering[clustering$Gene_Clusters==i,]$gene
  gene_entrez_id2GO <- clusterProfiler::bitr(gene, "SYMBOL", "ENTREZID", "org.Mm.eg.db", drop = TRUE)$ENTREZID    
  eG <- enrichGO(gene = gene_entrez_id2GO, #需要分析的基因的EntrezID
                 OrgDb = org.Mm.eg.db, #人基因数据库
                 pvalueCutoff =0.05, #设置pvalue界值
                 qvalueCutoff = 0.05, #设置qvalue界值(FDR校正后的p值）
                 ont="all", #选择功能富集的类型，可选BP、MF、CC，这里选择all。
                 readable =T)
  GO <- data.frame(eG@result)
  GO$group <- rep(i,nrow(GO))
  GO_new <- data.frame(Term=GO$Description,
                       qvalue=GO$qvalue,
                       pvalue=GO$pvalue,
                       Count=GO$Count,
                       Ontology=GO$ONTOLOGY,
                       gene=GO$geneID,
                       cluster=GO$group)
  GO_all <- rbind(GO_all,GO_new)
}
write.csv(GO_all,"SCT_Neu_BEAM_sig_gene_GO_select_cluster.csv")
GO <- read.csv("SCT_Neu_BEAM_sig_gene_GO_select_cluster.csv")

library(ggplot2)
library(forcats)

test <- GO[GO$cluster==1,]
p1<-ggplot(test)+
  geom_bar(aes(-log10(as.numeric(pvalue)),fct_reorder(Term,Ontology),fill = Ontology),stat="identity") + #将Term按照Ontology排序
  theme_bw()+
  theme(text=element_text(family='serif',size =10,face = 'bold'),
        axis.text = element_text(family='serif',size =10,face = 'bold',colour = 'black'),
        panel.grid = element_blank()) + 
  #geom_text(data=subset(test,qvalue<0.05,c('Count','Term')),aes(x=Count+3,y=as.factor(Term),label='*')) + #对于qvalue<0.05的term添加*
  labs(x='-log10(pvalue)',y='GO term',title = "SCT_Neu_cluster1_gene_GO_term") + 
  scale_fill_manual(values=c(BP = "#79B494", CC = "#D67E56", MF = "#848CBD"))
test <- GO[GO$cluster==2,]
p2<-ggplot(test)+
  geom_bar(aes(-log10(as.numeric(pvalue)),fct_reorder(Term,Ontology),fill = Ontology),stat="identity") + #将Term按照Ontology排序
  theme_bw()+
  theme(text=element_text(family='serif',size =10,face = 'bold'),
        axis.text = element_text(family='serif',size =10,face = 'bold',colour = 'black'),
        panel.grid = element_blank()) + 
  #geom_text(data=subset(test,qvalue<0.05,c('Count','Term')),aes(x=Count+3,y=as.factor(Term),label='*')) + #对于qvalue<0.05的term添加*
  labs(x='-log10(pvalue)',y='GO term',title = "SCT_Neu_cluster2_gene_GO_term") + 
  scale_fill_manual(values=c(BP = "#79B494", CC = "#D67E56", MF = "#848CBD"))
test <- GO[GO$cluster==3,]
p3<-ggplot(test)+
  geom_bar(aes(-log10(as.numeric(pvalue)),fct_reorder(Term,Ontology),fill = Ontology),stat="identity") + #将Term按照Ontology排序
  theme_bw()+
  theme(text=element_text(family='serif',size =10,face = 'bold'),
        axis.text = element_text(family='serif',size =10,face = 'bold',colour = 'black'),
        panel.grid = element_blank()) + 
  #geom_text(data=subset(test,qvalue<0.05,c('Count','Term')),aes(x=Count+3,y=as.factor(Term),label='*')) + #对于qvalue<0.05的term添加*
  labs(x='-log10(pvalue)',y='GO term',title = "SCT_Neu_cluster3_gene_GO_term") + 
  scale_fill_manual(values=c(BP = "#79B494", CC = "#D67E56", MF = "#848CBD"))
cowplot::plot_grid(p1,p2,p3,ncol = 1,align = 'hv')
ggsave("SCT_Neu_pseu_cluster_BEAM_GO_term_log10pvalue.pdf",width = 8,height = 8)

barplot(ego, title = paste0("GO terms enriched in ", i))

setwd("E:\\SCT\\data\\投稿\\plos\\Fix\\拟时序分析\\gene_pseu")
##单独的基因可视化
for (gene in BEAM_res_gene$gene_short_name){
  p=plot_cell_trajectory(cds, markers = gene ,use_color_gradient = T)+
    tidydr::theme_dr(xlength = 0.2, ylength = 0.2,arrow = arrow(length = unit(0.2, "inches"),type = "closed")) +
    theme(panel.grid = element_blank(),
          axis.title = element_text(face = 2,hjust = 0.03))
  ggsave(paste0(gene,"_Neu_pseu.pdf",sep=""),p,width = 6,height = 5)
}
dev.new()
for (gene in BEAM_res_gene$gene_short_name){
  p=plot_cell_trajectory(cds, markers = gene ,use_color_gradient = T)+
    tidydr::theme_dr(xlength = 0.2, ylength = 0.2,arrow = arrow(length = unit(0.2, "inches"),type = "closed")) +
    theme(panel.grid = element_blank(),
          axis.title = element_text(face = 2,hjust = 0.03))
  ggsave(paste0(gene,"_Neu_pseu.pdf",sep=""),p,width = 6,height = 5)
}

features = c("S100b","Ptgds","Snap25","Nefl","Spp1","Calb2")
pdf("SCT_Neu_pseu_cluster_gene_pseudotime.pdf",width = 6,height = 4)
plot_genes_in_pseudotime(cds[features,],
                         color_by = "group",
                         ncol = 3,
                         nrow = NULL)+scale_color_manual(values = c("#79B494", "#D67E56", "#848CBD"))
dev.off()
plot_genes_branched_pseudotime(cds[features,],
                               branch_point = 1,
                               color_by = "group",
                               method = 'loess',
                               cell_size=2,
                               ncol = 3)

##创建对象
# 获取 integrated assay 的数据和对应的元数据
integrated_data <- SeuratObject::LayerData(Neu, assay = "integrated")
integrated_cells <- colnames(integrated_data)
integrated_metadata <- Neu@meta.data[integrated_cells, ]

# 创建 monocle3 对象
cds <- monocle3::new_cell_data_set(
  expression_data = integrated_data,
  cell_metadata = integrated_metadata,
  gene_metadata = data.frame(
    gene_short_name = rownames(integrated_data), 
    row.names = rownames(integrated_data)
  )
)

cat("Successfully created monocle3 object with", ncol(cds), "cells and", nrow(cds), "genes\n")
# 数据预处理
cds <- preprocess_cds(cds, num_dim = 50)
# 批次校正
cds <- align_cds(cds, preprocess_method = "PCA")
# 降维
cds <- reduce_dimension(cds,reduction_method = "UMAP")
cds.embed <- cds@int_colData$reducedDims$UMAP
int.embed <- Embeddings(Neu, reduction = "umap")
int.embed <- int.embed[rownames(cds.embed), ]
cds@int_colData$reducedDims$UMAP <- int.embed
plot_cells(cds, label_groups_by_cluster = FALSE, color_cells_by = "group")
# 细胞聚类
cds <- cluster_cells(cds)
# 按分区可视化
plot_cells(cds, color_cells_by = "partition")
# 学习轨迹图
cds <- learn_graph(cds)
# 可视化轨迹图
plot_cells(cds, color_cells_by = "group", label_groups_by_cluster = FALSE, label_leaves = FALSE, label_branch_points = FALSE)
# 按伪时间排序细胞
cds <- order_cells(cds)
# 按伪时间可视化
plot_cells(cds, color_cells_by = "pseudotime")
genes <- monocle3::graph_test(cds, neighbor_graph = "principal_graph", reduction_method = "UMAP", cores = 32)
top10 <- genes %>%
  top_n(n = 10, morans_I) %>%
  pull(gene_short_name) %>%
  as.character()
top50 <- genes %>%
  top_n(n = 10, morans_I) %>%
  pull(gene_short_name) %>%
  as.character()
library(ClusterGVis)
library(scales)

plot_pseudotime_heatmap2(mat,
                         num_clusters = 3,
                         show_rownames = F,
                         return_heatmap = F,
                         hmcols = colorRampPalette(c("darkblue", "white","darkred"))(256))
# kmeans
# 首先检查 cds 对象类型
# 完全重新创建矩阵，避免任何继承的属性
mat <- cds@assays@data$counts
mat <- data.frame(mat)
mat <- mat[common_marker,]
mat <- mat[!grepl("^NA", rownames(mat)), ]
ck <- clusterData(mat,
                  cluster.method = "kmeans",
                  cluster.num =3)

# add line annotation
pdf('monocle3.pdf',height = 5,width = 6,onefile = F)
visCluster(object = ck,
           plot.type = "heatmap",
           ht.col.list = list(col_range=c(-2,0,2),
                              col_color=colorRampPalette(c("darkblue", "white","darkred"))(3)),
           #cluster.order = "group",
           add.sampleanno = F,
           show_row_names = F)
dev.off()
counts_matrix <- monocle3::exprs(cds)
mat <- data.frame(counts_matrix)
mat <- mat[common_marker,]
ck <- clusterData(mat,
                  cluster.method = "kmeans",
                  cluster.num = 3)
gene_select <- genes[genes$status=="OK",]
Idents(Neu) <- Neu$group
Idents(Neu) <- factor(Idents(Neu), levels = c("E17.5","P8","Adult"))

DotPlot(Neu,features = gene_select$gene_short_name[1:50])+RotatedAxis()+NoLegend()
ggsave("Neu_group_gene50.pdf",width = 20,height = 5)

avg_data <- AverageExpression(Neu,
                              features = gene_select$gene_short_name,
                              group.by = "group",
                     assays= "integrated")
mat_avg <- data.frame(avg_data)

# 寻找在E17.5高表达，随发育下降的基因
down_regulated <- mat_avg %>%
  filter(integrated.E17.5 > integrated.P8 & integrated.E17.5 > integrated.Adult&
           integrated.P8 > integrated.Adult) %>%
  arrange(desc(integrated.E17.5))

# 寻找在Adult高表达，随发育上升的基因  
up_regulated <- mat_avg %>%
  filter(integrated.Adult > integrated.P8 & integrated.Adult > integrated.E17.5 & 
           integrated.P8 > integrated.E17.5) %>%
  arrange(integrated.E17.5)



common_marker <- union(rownames(up_regulated),rownames(down_regulated))

library(COSG)
COSG_markers <- cosg(
  object1,
  groups = 'all',
  assay ='Spatial',
  slot='data',
  mu=1,
  n_genes_user=200)

celltype <- colnames(COSG_markers$names)
COSG_result_all <- data.frame(celltype=as.character(), markers=as.character(), COSG_score=as.numeric())

for (i in celltype) {
  COSG_result_new <- data.frame(
    celltype = rep(i, 200),
    markers = COSG_markers$names[[i]][1:200],  # 使用 [[i]] 而不是 $i
    COSG_score = COSG_markers$scores[[i]][1:200]  # 同样使用 [[i]]
  )
  COSG_result_all <- rbind(COSG_result_all, COSG_result_new)
}

Neu_marker_E <- COSG_result_all[COSG_result_all$celltype=="SGN",]$markers

COSG_markers <- cosg(
  object2,
  groups = 'all',
  assay ='Spatial',
  slot='data',
  mu=1,
  n_genes_user=200)

celltype <- colnames(COSG_markers$names)
COSG_result_all <- data.frame(celltype=as.character(), markers=as.character(), COSG_score=as.numeric())

for (i in celltype) {
  COSG_result_new <- data.frame(
    celltype = rep(i, 200),
    markers = COSG_markers$names[[i]][1:200],  # 使用 [[i]] 而不是 $i
    COSG_score = COSG_markers$scores[[i]][1:200]  # 同样使用 [[i]]
  )
  COSG_result_all <- rbind(COSG_result_all, COSG_result_new)
}

Neu_marker_P8 <- COSG_result_all[COSG_result_all$celltype=="SGN",]$markers

COSG_markers <- cosg(
  object3,
  groups = 'all',
  assay ='Spatial',
  slot='data',
  mu=1,
  n_genes_user=200)

celltype <- colnames(COSG_markers$names)
COSG_result_all <- data.frame(celltype=as.character(), markers=as.character(), COSG_score=as.numeric())

for (i in celltype) {
  COSG_result_new <- data.frame(
    celltype = rep(i, 200),
    markers = COSG_markers$names[[i]][1:200],  # 使用 [[i]] 而不是 $i
    COSG_score = COSG_markers$scores[[i]][1:200]  # 同样使用 [[i]]
  )
  COSG_result_all <- rbind(COSG_result_all, COSG_result_new)
}

Neu_marker_Adult <- COSG_result_all[COSG_result_all$celltype=="SGN",]$markers

Neu_marker <- union(Neu_marker_E,Neu_marker_P8)
Neu_marker <- union(Neu_marker,Neu_marker_Adult)
