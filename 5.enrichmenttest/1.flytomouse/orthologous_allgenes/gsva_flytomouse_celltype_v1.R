callback = function(hc, mat){
  sv = svd(t(mat))$v[,1]
  dend = reorder(as.dendrogram(hc), wts = sv)
  as.hclust(dend)
}

library(pheatmap)
library(plyr)
library(dplyr)

flyObj.markers <- readRDS('../../celltype_v1.Rds')
flyObj.markers <- subset(flyObj.markers, p_val_adj <= 0.05)
head(flyObj.markers); nrow(flyObj.markers) # 2703

counts <- data.frame(table(flyObj.markers$gene))
reduns <- subset(counts, Freq > 1); head(reduns); nrow(reduns) # 494
flyObj.markers <- subset(flyObj.markers, !gene %in% as.character(reduns$Var1))
head(flyObj.markers); nrow(flyObj.markers) # 1572


conserved <- read.delim('../../../2.genelists/fly/diopt_flytomouse.usegenes.txt')
conserved <- conserved[, c(1,8)]
rownames(conserved) <- conserved[, 1]
head(conserved); nrow(conserved) # 6470
flyObj.markers <- subset(flyObj.markers, gene %in% conserved$FlyGeneID); nrow(flyObj.markers) # 1155
summary(flyObj.markers$cluster)
# PH 1    PH    PM    PM 120    LM    CC    GST-rich    Adipohemocyte  
# 464     41    32    47        333   66    78          94 
#write.table(flyObj.markers, 'gsva_flytomouse_celltype_v1.orthologous.txt', row.names = F, col.names = T, sep = '\t', quote = F)


### Gene set ###
library(GSEABase)
library(GSVA)

gsInfo <- data.frame(row.names = levels(flyObj.markers$cluster), Source = levels(flyObj.markers$cluster)); head(gsInfo)
gsInfo <- AnnotatedDataFrame(data=gsInfo); gsInfo

targetgene_w <- data.frame(matrix(nrow = max(summary(flyObj.markers$cluster)), ncol = length(levels(flyObj.markers$cluster))))
colnames(targetgene_w) <- levels(flyObj.markers$cluster)

gs_list <- list()
head(conserved, n = 3)
for (i in c(1:length(levels(flyObj.markers$cluster)))) {
  celltype <- levels(flyObj.markers$cluster)[i]
  refgene <- as.character(subset(flyObj.markers, cluster == celltype)$gene)
  #targetgene <- as.character(subset(conserved, FlyGeneID %in% refgene)$MouseSymbol)
  targetgene <- as.character(conserved[refgene, 'MouseSymbol'])
  if (length(targetgene) > max(summary(flyObj.markers$cluster))) {
    targetgene <- targetgene[1:max(summary(flyObj.markers$cluster))]
  }
  print (length(targetgene))
  gs_list[[i]] <- GeneSet(targetgene, setName = celltype)
  if (length(targetgene) < max(summary(flyObj.markers$cluster))) {
    targetgene <- c(targetgene, rep('', times = max(summary(flyObj.markers$cluster))-length(targetgene)))
  }
  targetgene_w[, celltype] <- targetgene
}
head(targetgene_w)
#write.table(targetgene_w, 'gsva_flytomouse_celltype_v1.orthologous.allgenes.MouseIDs.txt', row.names = F, col.names = T, sep = '\t', quote = F)
mySignatures <- GeneSetCollection(object = gs_list) #


### GSVA ###
# MCA
mouseexpr <- readRDS('../../../mouse/mca/tmp/mca_pseudobulk_noNeu.Rds')
mouseexpr <- mouseexpr[intersect(as.character(conserved$MouseSymbol), rownames(mouseexpr)), ]; dim(mouseexpr) # 5470   10
mouseexpr$sums <- rowSums(mouseexpr)
mouseexpr <- subset(mouseexpr, sums > 0 )
mouseexpr$sums <- NULL
head(mouseexpr); dim(mouseexpr) # 5463   10

label <- data.frame(row.names = colnames(mouseexpr), celltype = colnames(mouseexpr)); label
label <- AnnotatedDataFrame(data=label)

exprSet <- ExpressionSet(as.matrix(mouseexpr), phenoData = label, annotation = 'Symbol'); exprSet

es.gsva <- gsva(expr = exprSet, gset.idx.list = mySignatures,
                method='gsva', verbose=T, kcdf = 'Gaussian', parallel.sz =1)
gsvaScores <- as.data.frame(t(exprs(es.gsva)))
head(gsvaScores)
gsvaScores_w <- data.frame(celltype = rownames(gsvaScores), gsvaScores, check.rows = F, check.names = F)
#write.table(gsvaScores_w, 'gsva_flytomouse_celltype.mca.txt', quote = F, sep = '\t', row.names = F, col.names = T)
min(gsvaScores); max(gsvaScores)

pheatmap(gsvaScores, cluster_rows = F, cluster_cols = F, scale = 'none', treeheight_row = 10, treeheight_col = 10,
         cellwidth = 10, cellheight = 10, fontsize_row = 9, fontsize_col = 9, border_color = NA,
         clustering_distance_rows = 'euclidean', clustering_distance_cols = 'euclidean', clustering_method = 'complete', clustering_callback = callback, 
         breaks = c(-0.6, -0.5, -0.4, -0.3, -0.2, -0.1, 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6),
         color = colorRampPalette(c("#3a5fcd", 'white', "#ee0000"))(12), filename = 'gsva_flytomouse_celltype.mca.pdf')


# TM 10X
mouseexpr <- readRDS('../../../mouse/tabulamuris_10x/tmp/tm_10x_pseudobulk.Rds')
mouseexpr <- mouseexpr[intersect(as.character(conserved$MouseSymbol), rownames(mouseexpr)), ]; dim(mouseexpr) # 5493   11
mouseexpr$sums <- rowSums(mouseexpr)
mouseexpr <- subset(mouseexpr, sums > 0 )
mouseexpr$sums <- NULL
head(mouseexpr); dim(mouseexpr) # 5493   11

label <- data.frame(row.names = colnames(mouseexpr), celltype = colnames(mouseexpr)); label
label <- AnnotatedDataFrame(data=label)

exprSet <- ExpressionSet(as.matrix(mouseexpr), phenoData = label, annotation = 'Symbol'); exprSet

es.gsva <- gsva(expr = exprSet, gset.idx.list = mySignatures,
                method='gsva', verbose=T, kcdf = 'Gaussian', parallel.sz =1)
gsvaScores <- as.data.frame(t(exprs(es.gsva)))
head(gsvaScores)
gsvaScores_w <- data.frame(celltype = rownames(gsvaScores), gsvaScores, check.rows = F, check.names = F)
#write.table(gsvaScores_w, 'gsva_flytomouse_celltype.tm10x.txt', quote = F, sep = '\t', row.names = F, col.names = T)
min(gsvaScores); max(gsvaScores)

pheatmap(gsvaScores, cluster_rows = F, cluster_cols = F, scale = 'none', treeheight_row = 10, treeheight_col = 10,
         cellwidth = 10, cellheight = 10, fontsize_row = 9, fontsize_col = 9, border_color = NA,
         clustering_distance_rows = 'euclidean', clustering_distance_cols = 'euclidean', clustering_method = 'complete', clustering_callback = callback, 
         breaks = c(-0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1, 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7),
         color = colorRampPalette(c("#3a5fcd", 'white', "#ee0000"))(14), filename = 'gsva_flytomouse_celltype.tm10x.pdf')


# TM SS2
mouseexpr <- readRDS('../../../mouse/tabulamuris_ss2/tmp/tm_ss2_pseudobulk.Rds')
mouseexpr <- mouseexpr[intersect(as.character(conserved$MouseSymbol), rownames(mouseexpr)), ]; dim(mouseexpr) # 6027   11
mouseexpr$sums <- rowSums(mouseexpr)
mouseexpr <- subset(mouseexpr, sums > 0 )
mouseexpr$sums <- NULL
head(mouseexpr); dim(mouseexpr) # 6027   11

label <- data.frame(row.names = colnames(mouseexpr), celltype = colnames(mouseexpr)); label
label <- AnnotatedDataFrame(data=label)

exprSet <- ExpressionSet(as.matrix(mouseexpr), phenoData = label, annotation = 'Symbol'); exprSet

es.gsva <- gsva(expr = exprSet, gset.idx.list = mySignatures,
                method='gsva', verbose=T, kcdf = 'Gaussian', parallel.sz =1)
gsvaScores <- as.data.frame(t(exprs(es.gsva)))
head(gsvaScores)
gsvaScores_w <- data.frame(celltype = rownames(gsvaScores), gsvaScores, check.rows = F, check.names = F)
#write.table(gsvaScores_w, 'gsva_flytomouse_celltype.tmss2.txt', quote = F, sep = '\t', row.names = F, col.names = T)
min(gsvaScores); max(gsvaScores)

pheatmap(gsvaScores, cluster_rows = F, cluster_cols = F, scale = 'none', treeheight_row = 10, treeheight_col = 10,
         cellwidth = 10, cellheight = 10, fontsize_row = 9, fontsize_col = 9, border_color = NA,
         clustering_distance_rows = 'euclidean', clustering_distance_cols = 'euclidean', clustering_method = 'complete', clustering_callback = callback, 
         breaks = c(-0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1, 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7),
         color = colorRampPalette(c("#3a5fcd", 'white', "#ee0000"))(14), filename = 'gsva_flytomouse_celltype.tmss2.pdf')
dev.off()


