######## correlation analysis ###########
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
##新注释空间数据集调整
object1_anno <- readRDS("E17_SCT_R_anno.Rds")
levels(object1_anno)
ids <- c("OS","MaC","NeC","EryC","Ube2c+",
         "SLg_FC","SLb_FC",
         "SGN",
         "SchC",
         "RMC",
         "Stria",
         "OSC",
         "KO",
         "HC",
         "DC/PC")
names(ids) <- levels(object1_anno)
object1_anno <- RenameIdents(object1_anno,ids)
object1_anno$cellType <- Idents(object1_anno)
COSG_markers <- cosg(
  object1_anno,
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
saveRDS(object1_anno,"E17_SCT_R_RScore.Rds")
write.csv(COSG_result_all,"E17_SCT_R_marker_COSG_score.csv")

object3_anno <- readRDS("P8_SCT_Up_anno.Rds")
levels(object3_anno)
ids <- c("OS","NeC","OsC",
         "SLg_FC","SLb_FC",
         "SGN",
         "SchC",
         "RMC",
         "Stria","SP",
         "OSC",
         "IdC","ISC","IBC","IPhC",
         "TBC",
         "HC",
         "DC/PC")
names(ids) <- levels(object3_anno)
object3_anno <- RenameIdents(object3_anno,ids)
object3_anno$cellType <- Idents(object3_anno)
COSG_markers <- cosg(
  object3_anno,
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
saveRDS(object3_anno,"P8_SCT_Up_RScore.Rds")
write.csv(COSG_result_all,"P8_SCT_Up_marker_COSG_score.csv")

##计算相关性函数
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
  
  geTable = base::merge(Sp1,Sp2, by='row.names', all=F)
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

save_pheatmap_pdf <- function(x, filename, width=7, height=7) {
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  pdf(filename, width=width, height=height)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}
###
scRNA <- readRDS("E16_scRNA_RScore.Rds")
spRNA <- readRDS("E17_SCT_R_RScore.Rds")
ref_gene <- read.csv("E16_marker_COSG_score.csv")
que_gene <- read.csv("E17_SCT_R_marker_COSG_score.csv")
ref <- list(c(scRNA,"RNA","cellType","scRNA"))
que <- list(c(spRNA,"Spatial","cellType","spatial"))
feature <- union(ref_gene$markers,que_gene$markers)
tmp <- cortable(ref,que,features = feature)
Corr.Coeff.Table <- tmp[[1]]
p.value.table <- tmp[[2]]
write.csv(data.frame(Corr.Coeff.Table),"E16_scRNA_E17_R_spatial_coeff_corration.csv")
write.csv(data.frame(p.value.table),"E16_scRNA_E17_R_spatial_pvalue_corration.csv")

scRNA <- readRDS("P7_scRNA_RScore.Rds")
spRNA <- readRDS("P8_SCT_Up_RScore.Rds")
ref_gene <- read.csv("P7_marker_COSG_score.csv")
que_gene <- read.csv("P8_SCT_Up_marker_COSG_score.csv")
ref <- list(c(scRNA,"RNA","cellType","scRNA"))
que <- list(c(spRNA,"Spatial","cellType","spatial"))
feature <- union(ref_gene$markers,que_gene$markers)
tmp <- cortable(ref,que,features = feature)
Corr.Coeff.Table <- tmp[[1]]
p.value.table <- tmp[[2]]
write.csv(data.frame(Corr.Coeff.Table),"P7_scRNA_P8_Up_spatial_coeff_corration.csv")
write.csv(data.frame(p.value.table),"P7_scRNA_P8_Up_spatial_pvalue_corration.csv")

cor_matrix <- read.csv("E16_scRNA_E17_R_spatial_coeff_corration.csv")
rownames(cor_matrix) <- cor_matrix$X
cor_matrix<- cor_matrix[,-1]
cor_matrix <- as.matrix(cor_matrix)
p <- pheatmap(t(cor_matrix),scale="none",show_colnames = T,show_rownames = T,fontsize = 10,
              cluster_rows = F,cluster_cols = F,
              border=F,border_color = "NA",
              color = colorRampPalette(c("darkblue", "white","darkred"))(256),
              cellwidth = 20,cellheight = 16,
              #gaps_row = c(1,2,3,8,11,12),
              main = "E16_scRNA VS E17_L_spatial",
              fontsize_number = 12, number_color = "black")
save_pheatmap_pdf(p,"E16_scRNA VS E17_L_spatial Pearson correlation coefficient.pdf",8,6)

cor_matrix <- read.csv("P7_scRNA_P8_Up_spatial_coeff_corration.csv")
rownames(cor_matrix) <- cor_matrix$X
cor_matrix<- cor_matrix[,-1]
cor_matrix <- as.matrix(cor_matrix)
p<-pheatmap(t(cor_matrix),scale="none",show_colnames = T,show_rownames = T,fontsize = 10,
            cluster_rows = F,cluster_cols = F,
            border=F,border_color = "NA",
            color = colorRampPalette(c("darkblue", "white","darkred"))(256),
            cellwidth = 20,cellheight = 16,
            #gaps_row = c(1,2,3,8,11,12),
            main = "P7_scRNA VS P8_Up_spatial",
            fontsize_number = 12, number_color = "black")

save_pheatmap_pdf(p, "P7_scRNA VS P8_Up_spatial Pearson correlation coefficient_top50.pdf",10,6)


##coreelation analysis beween biological replicates
spRNA1 <- readRDS("E17_spatial_RNA_RScore.Rds")
spRNA2 <- readRDS("E17_SCT_R_RScore.Rds")
ref_gene <- read.csv("E17_spatial_marker_COSG_score.csv")
que_gene <- read.csv("Fig1\\E17_SCT_R_marker_COSG_score.csv")
ref <- list(c(spRNA1,"Spatial","cellType","spatial"))
que <- list(c(spRNA2,"Spatial","cellType","spatial"))
feature <- union(ref_gene$markers,que_gene$markers)
tmp <- cortable(ref,que,features = feature)
Corr.Coeff.Table <- tmp[[1]]
p.value.table <- tmp[[2]]
write.csv(data.frame(Corr.Coeff.Table),"E17_L_spatial_E17_R_spatial_coeff_corration.csv")
write.csv(data.frame(p.value.table),"E17_L_spatial_E17_R_spatial_pvalue_corration.csv")

spRNA1 <- readRDS("P8_spatial_RNA_RScore.Rds")
spRNA2 <- readRDS("P8_SCT_Up_RScore.Rds")
ref_gene <- read.csv("P8_spatial_marker_COSG_score.csv")
que_gene <- read.csv("P8_SCT_Up_marker_COSG_score.csv")
ref <- list(c(spRNA1,"Spatial","cellType","spatial"))
que <- list(c(spRNA2,"Spatial","cellType","spatial"))
feature <- union(ref_gene$markers,que_gene$markers)
tmp <- cortable(ref,que,features = feature)
Corr.Coeff.Table <- tmp[[1]]
p.value.table <- tmp[[2]]
write.csv(data.frame(Corr.Coeff.Table),"P8_down_P8_Up_spatial_coeff_corration.csv")
write.csv(data.frame(p.value.table),"P8_down_P8_Up_spatial_pvalue_corration.csv")

cor_matrix <- read.csv("E17_L_spatial_E17_R_spatial_coeff_corration.csv")
rownames(cor_matrix) <- cor_matrix$X
cor_matrix<- cor_matrix[,-1]
cor_matrix <- as.matrix(cor_matrix)
p <- pheatmap(t(cor_matrix),scale="none",show_colnames = T,show_rownames = T,fontsize = 10,
              cluster_rows = F,cluster_cols = F,
              border=F,border_color = "NA",
              color = colorRampPalette(c("darkblue", "white","darkred"))(256),
              cellwidth = 20,cellheight = 16,
              #gaps_row = c(1,2,3,8,11,12),
              main = "E17_R_scRNA VS E17_L_spatial",
              fontsize_number = 12, number_color = "black")
save_pheatmap_pdf(p,"E17_R VS E17_L_spatial Pearson correlation coefficient.pdf",8,6)

cor_matrix <- read.csv("P8_down_P8_Up_spatial_coeff_corration.csv")
rownames(cor_matrix) <- cor_matrix$X
cor_matrix<- cor_matrix[,-1]
cor_matrix <- as.matrix(cor_matrix)
p<-pheatmap(t(cor_matrix),scale="none",show_colnames = T,show_rownames = T,fontsize = 10,
            cluster_rows = F,cluster_cols = F,
            border=F,border_color = "NA",
            color = colorRampPalette(c("darkblue", "white","darkred"))(256),
            cellwidth = 20,cellheight = 16,
            #gaps_row = c(1,2,3,8,11,12),
            main = "P8_Down VS P8_Up_spatial",
            fontsize_number = 12, number_color = "black")

save_pheatmap_pdf(p, "P8_Down VS P8_Up_spatial Pearson correlation coefficient_top50.pdf",10,6)
