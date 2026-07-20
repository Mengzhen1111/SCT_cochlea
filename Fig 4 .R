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
library(Cairo)
##载入数据
object1 <- readRDS("E17_SCT_L_anno.Rds")
object2 <- readRDS("P8_SCT_Down_anno.Rds")
object3 <- readRDS("Adult_SCT_Right_anno.Rds")
object1$group <- "E17.5"
object2$group <- "P8"
object3$group <- "Adult"

data.list <- list(E17=object1,P8=object2,Adult=object3)
for (i in 1:length(data.list)) {
  data.list[[i]] <- NormalizeData(data.list[[i]], normalization.method = "LogNormalize", scale.factor = 10000)
  data.list[[i]] <- FindVariableFeatures(data.list[[i]], selection.method = "vst", nfeatures = 2000)}
features <- SelectIntegrationFeatures(object.list = data.list)
data.anchors <- FindIntegrationAnchors(object.list = data.list, anchor.features = features)
data_in <- IntegrateData(anchorset = data.anchors)
DefaultAssay(data_in) <- "integrated"
data_in <- ScaleData(data_in, verbose = FALSE)
data_in <- RunPCA(data_in, npcs = 100, verbose = FALSE)
data_in <- RunUMAP(data_in, reduction = "pca", dims = 1:50)
data_in <- FindNeighbors(data_in, reduction = "pca", dims = 1:50)
data_in <- FindClusters(data_in, resolution = 0.8,verbose = FALSE)
DimPlot(data_in, reduction = "umap", group.by = c("ident", "group"))
Idents(data_in) <- data_in$cellType
Neu <- subset(data_in,idents = "SGN")
HC <- subset(data_in,idents = c("HC","OHC","IHC"))

saveRDS(Neu,"SCT_group_neu.Rds")
saveRDS(HC,"SCT_group_HC.Rds")

HC <- readRDS("SCT_group_HC.Rds")
Idents(HC)<-HC$group

HC <- ScaleData(HC, verbose = FALSE)
HC <- RunPCA(HC, npcs = 50, verbose = FALSE)
HC <- RunUMAP(HC, reduction = "pca", dims = 1:50)
HC <- FindNeighbors(HC, reduction = "pca", dims = 1:50)
HC <- FindClusters(HC, resolution = 0.1,verbose = FALSE)
DimPlot(HC, reduction = "umap", group.by = c("group"))
Idents(HC)<-HC$group

Idents(HC) <- factor(Idents(HC), levels = c("E17.5","P8","Adult"))

DimPlot(HC, reduction = "umap",pt.size = 2, cols = c("#FBC04C","#52C7B7", "#4C77C6"),label = F,label.size = 5,group.by = c("group"))+
  tidydr::theme_dr(xlength = 0.2, ylength = 0.2,arrow = arrow(length = unit(0.2, "inches"),type = "closed")) +
  theme(panel.grid = element_blank(),
        axis.title = element_text(face = 2,hjust = 0.03))
ggsave("SCT_HC_group_Umap.pdf",width = 5.5,height = 5)


object1 <- readRDS("E17_SCT_L_anno.Rds")
object2 <- readRDS("P8_SCT_Down_anno.Rds")
object3 <- readRDS("Adult_SCT_Right_anno.Rds")
E17 <- subset(object1,idents="HC")
E17$group <- "E17.5"
P8 <- subset(object2,idents=c("OHC","IHC"))
P8$group <- "P8"
Adult <- subset(object3,idents=c("OHC","IHC"))
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
ordering_genes <- rownames(subset(diff_test_res, pval < 0.05))

cds <- setOrderingFilter(cds,ordering_genes)
plot_ordering_genes(cds)

##降维

cds <- reduceDimension(cds,method='DDRTree')
##排序
cds <- orderCells(cds)
# 可视化轨迹，按伪时间着色
plot_cell_trajectory(cds, color_by = "Pseudotime")+
  tidydr::theme_dr(xlength = 0.2, ylength = 0.2,arrow = arrow(length = unit(0.2, "inches"),type = "closed")) +
  theme(panel.grid = element_blank(),
        axis.title = element_text(face = 2,hjust = 0.03))
ggsave("SCT_HC_pseu_Pseudotime.pdf",width = 5,height = 3)

# 可视化轨迹，按样本着色
plot_cell_trajectory(cds, color_by = "group")+scale_colour_manual(values = c("#79B494", "#D67E56",  "#848CBD"))+
  tidydr::theme_dr(xlength = 0.2, ylength = 0.2,arrow = arrow(length = unit(0.2, "inches"),type = "closed")) +
  theme(panel.grid = element_blank(),
        axis.title = element_text(face = 2,hjust = 0.03))
ggsave("SCT_HC_pseu_group.pdf",width = 5,height = 3)
# 可视化轨迹，按样本着色
plot_cell_trajectory(cds, color_by = "State")+
  tidydr::theme_dr(xlength = 0.2, ylength = 0.2,arrow = arrow(length = unit(0.2, "inches"),type = "closed")) +
  theme(panel.grid = element_blank(),
        axis.title = element_text(face = 2,hjust = 0.03))
ggsave("SCT_HC_pseu_State.pdf",width = 5,height = 3)

##拟时序差异基因热图绘制
cds_DGT_pseudotimegenes <- differentialGeneTest(cds, fullModelFormulaStr = "~group")
cds_DGT_pseudotimegenes_sig <- subset(cds_DGT_pseudotimegenes, pval < 0.05)
write.csv(cds_DGT_pseudotimegenes_sig,"SCT_HC_peudotimegenes_pvalue_0.05.csv")
cds_DGT_pseudotimegenes_sig_min <- subset(cds_DGT_pseudotimegenes, pval < 0.05&num_cells_expressed>30)

# 然后绘制热图
plot_pseudotime_heatmap(cds[cds_DGT_pseudotimegenes_sig_min$gene_short_name,], 
                        num_cluster = 3, 
                        show_rownames = F, 
                        return_heatmap = T,
                        #add_annotation_col = annotation_col,
                        hmcols = colorRampPalette(c("navy","white","firebrick3"))(100))


plot_pseudotime_heatmap(cds[Time_genes,], 
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

p=plot_pseudotime_heatmap(cds[cds_DGT_pseudotimegenes_sig_min$gene_short_name,], 
                          num_cluster = 3, 
                          show_rownames = T, 
                          return_heatmap = T,
                          hmcols = colorRampPalette(c("navy","white","firebrick3"))(100))

clusters <- cutree(p$tree_row, k = 3)
clustering <- data.frame(clusters)
clustering[,1] <- as.character(clustering[,1])
colnames(clustering) <- "Gene_Clusters"
clustering$gene <- rownames(clustering)
##分别进行GO分析

GO_all<-data.frame(Description=character(),
                   qvalue=numeric(),
                   pvalue=numeric(),
                   RichFactor=numeric(),
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
                       RichFactor=GO$RichFactor,
                       Count=GO$Count,
                       Ontology=GO$ONTOLOGY,
                       gene=GO$geneID,
                       cluster=GO$group)
  GO_all <- rbind(GO_all,GO_new)
}
write.csv(GO_all,"SCT_HC_pseu_sig_gene_GO_select_cluster.csv")

HC_pseu_GO_gene <- unique(unlist(strsplit(GO$gene, "/")))

HC_pseu_GO_gene_sig <- cds_DGT_pseudotimegenes_sig[cds_DGT_pseudotimegenes_sig$gene_short_name%in%HC_pseu_GO_gene,]
write.csv(HC_pseu_GO_gene_sig,"HC_pseu_GO_gene_sig.csv")

pdf("SCT_HC_pseudotime_heatmap.pdf",height = 8,width = 7)
plot_pseudotime_heatmap(cds[HC_pseu_GO_gene,], 
                          num_cluster = 3, 
                          show_rownames = T, 
                          return_heatmap = T,
                          hmcols = colorRampPalette(c("navy","white","firebrick3"))(100))
dev.off()
dev.new()
library(ggplot2)
library(forcats)
GO_last <- read.csv("SCT_HC_pseu_sig_gene_GO_select_cluster_results_heatmap.csv")
GO <- GO_all[GO_all$Term %in% GO_last$Term, ]
##Style1
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
ggsave("SCT_HC_pseu_cluster_GO_term_genecount.pdf",width = 8,height = 15)

##Style2
mytheme <- theme(axis.title=element_text(face="bold", size=14,colour = 'black'), #坐标轴标题
                 axis.text=element_text(face="bold", size=14,colour = 'black'), #坐标轴标签
                 axis.line = element_line(size=0.5, colour = 'black'), #轴线
                 panel.background = element_rect(color='black'), #绘图区边框
                 legend.key = element_blank() #关闭图例边框
)

#绘制GO气泡图
test <- GO[GO$cluster==1,]
p1 <-ggplot(test,aes(x=-1*log10(as.numeric(pvalue)),y=Term,colour=RichFactor,size=Count))+
  geom_point()+
  scale_size(range=c(2, 8))+
  scale_colour_gradient(low = "blue",high = "red")+
  theme_bw()+
  ylab("GO Terms")+
  xlab("-1*log10(pvalue)")+labs(title = "Cluster1")+
  labs(color=expression(-log[10](PValue)))+
  theme(legend.title=element_text(size=14), legend.text = element_text(size=14))+
  theme(axis.title.y = element_text(margin = margin(r = 50)),axis.title.x = element_text(margin = margin(t = 20)))+
  theme(axis.text.x = element_text(face ="bold",color="black",angle=0,vjust=1))+mytheme

test <- GO[GO$cluster==2,]
p2 <-ggplot(test,aes(x=-1*log10(as.numeric(pvalue)),y=Term,colour=RichFactor,size=Count))+
  geom_point()+
  scale_size(range=c(2, 8))+
  scale_colour_gradient(low = "blue",high = "red")+
  theme_bw()+
  ylab("GO Terms")+
  xlab("-1*log10(pvalue)")+labs(title = "Cluster2")+
  labs(color=expression(-log[10](PValue)))+
  theme(legend.title=element_text(size=14), legend.text = element_text(size=14))+
  theme(axis.title.y = element_text(margin = margin(r = 50)),axis.title.x = element_text(margin = margin(t = 20)))+
  theme(axis.text.x = element_text(face ="bold",color="black",angle=0,vjust=1))+mytheme
test <- GO[GO$cluster==3,]
p3 <-ggplot(test,aes(x=-1*log10(as.numeric(pvalue)),y=Term,colour=RichFactor,size=Count))+
  geom_point()+
  scale_size(range=c(2, 8))+
  scale_colour_gradient(low = "blue",high = "red")+
  theme_bw()+
  ylab("GO Terms")+
  xlab("-1*log10(pvalue)")+labs(title = "Cluster3")+
  labs(color=expression(-log[10](PValue)))+
  theme(legend.title=element_text(size=14), legend.text = element_text(size=14))+
  theme(axis.title.y = element_text(margin = margin(r = 50)),axis.title.x = element_text(margin = margin(t = 20)))+
  theme(axis.text.x = element_text(face ="bold",color="black",angle=0,vjust=1))+mytheme
cowplot::plot_grid(p1,p2,p3,ncol = 1,align = 'hv')

ggsave("SCT_HC_pseu_cluster_GO_term_barplot_cluster1.pdf",p1,width = 12,height =8)
ggsave("SCT_HC_pseu_cluster_GO_term_barplot_cluster2.pdf",p2,width = 10,height =4)
ggsave("SCT_HC_pseu_cluster_GO_term_barplot_cluster3.pdf",p3,width = 13,height =4)


#Monocle2提供了一种特殊的统计检验方法：branched expression analysis modeling（BEAM），可以对不同的分支事件进行分析。

saveRDS(cds, file = "SCT_HC_pseudotime_analysis_results.Rds")

##分支节点基因热图
BEAM_res <- BEAM(cds, branch_point = 1, cores = 1,progenitor_method = "duplicate")
BEAM_res <- BEAM_res[order(BEAM_res$qval),]
BEAM_res <- BEAM_res[,c("gene_short_name", "pval", "qval")]
BEAM_res_gene <- subset(BEAM_res,pval < 0.05)

write.csv(BEAM_res_gene,"SCT_HC_pseudotime_BEAM_sig_gene.csv")
plot_genes_branched_heatmap(cds[BEAM_res_gene$gene_short_name,],
                            branch_point = 1,
                            num_clusters = 3,
                            cores = 1,
                            use_gene_short_name = T,
                            show_rownames = T)

pdf("SCT_HC_pseudotime_branched_heatmap.pdf",height = 8,width = 6)
pdf("SCT_HC_branched_heatmap.pdf",height = 8,width = 6)

plot_genes_branched_heatmap(cds[BEAM_res_gene$gene_short_name,],
                            branch_point = 1,
                            num_clusters = 3,
                            cores = 1,
                            show_rownames = T,
                            #return_heatmap = T,
                            use_gene_short_name = T)
dev.off()
dev.new()
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
                   RichFator=numeric(),
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
                       RichFactor=GO$RichFactor,
                       Count=GO$Count,
                       Ontology=GO$ONTOLOGY,
                       gene=GO$geneID,
                       cluster=GO$group)
  GO_all <- rbind(GO_all,GO_new)
}
write.csv(GO_all,"SCT_HC_BEAM_sig_gene_GO_select_cluster.csv")
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

##单独的基因可视化
for (gene in BEAM_res_gene$gene_short_name){
  p=plot_cell_trajectory(cds, markers = gene ,use_color_gradient = T,cell_size = 2)+
    tidydr::theme_dr(xlength = 0.2, ylength = 0.2,arrow = arrow(length = unit(0.2, "inches"),type = "closed")) +
    theme(panel.grid = element_blank(),
          axis.title = element_text(face = 2,hjust = 0.03))
  ggsave(paste0(gene,"_HC_pseu.pdf",sep=""),p,width = 6,height = 4)
}


features = c("Rpl26","Tgfb2","Tuba1a","Rpl4",
             "Npy","Tectb","Otog","Calb1",
             "Calb2","Ocm","Gsn","Calm2")
pdf("SCT_HC_pseu_cluster_gene_pseudotime.pdf",width = 9,height = 6)
plot_genes_in_pseudotime(cds[features,],
                         color_by = "group",
                         panel_order =c("Rpl26","Tgfb2","Tuba1a","Rpl4",
                                        "Npy","Tectb","Otog","Calb1",
                                        "Calb2","Ocm","Gsn","Calm2"),
                         ncol = 4,
                         cell_size = 2,
                         nrow = NULL)+scale_color_manual(values = c("#79B494", "#D67E56", "#848CBD"))
dev.off()
plot_genes_branched_pseudotime(cds[features,],
                               branch_point = 1,
                               color_by = "group",
                               panel_order =c("Rpl26","Tgfb2","Tuba1a","Rpl4",
                                              "Npy","Tectb","Otog","Calb1",
                                              "Calb2","Ocm","Gsn","Calm2"),
                               method = 'loess',
                               cell_size=2,
                               ncol = 4)


DefaultAssay(HC)<-"Spatial"
Idents(HC) <- HC$group
DotPlot(HC, features = c("Rpl26","Tgfb2","Tuba1a","Rpl4",
                          "Npy","Tectb","Otog","Calb1",
                          "Calb2","Ocm","Gsn","Calm2"))+
  RotatedAxis()+coord_flip()
ggsave("SCT_HC_pseu_gene_dotplot.pdf",width = 4,height = 6)
