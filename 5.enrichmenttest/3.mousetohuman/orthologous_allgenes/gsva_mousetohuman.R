callback = function(hc, mat){
  sv = svd(t(mat))$v[,1]
  dend = reorder(as.dendrogram(hc), wts = sv)
  as.hclust(dend)
}

library(pheatmap)
library(plyr)
library(dplyr)

mca.markers <- read.delim('../../../mouse/mca/degs/markers.mca.celltype_refined_noNeu.MAST.txt')
mca.markers <- subset(mca.markers, p_val_adj <= 0.05)
mca.markers$cluster <- factor(mca.markers$cluster, levels = unique(mca.markers$cluster))
mca.markers$cluster <- mapvalues(mca.markers$cluster, from = levels(mca.markers$cluster), to = paste(levels(mca.markers$cluster), '(MCA)', sep = ' '))
genecounts <- data.frame(table(mca.markers$gene))
nonreduns <- subset(genecounts, Freq == 1); head(nonreduns); nrow(nonreduns) # 723
mca.markers <- subset(mca.markers, gene %in% nonreduns$Var1)
mca.markers <- droplevels(mca.markers); nrow(mca.markers) # 723

tm10x.markers <- read.delim('../../../mouse/tabulamuris_10x/degs/markers.mouse.10x.celltype_refined.MAST.txt')
tm10x.markers <- subset(tm10x.markers, p_val_adj <= 0.05)
tm10x.markers$cluster <- factor(tm10x.markers$cluster, levels = unique(tm10x.markers$cluster))
tm10x.markers$cluster <- mapvalues(tm10x.markers$cluster, from = levels(tm10x.markers$cluster), to = paste(levels(tm10x.markers$cluster), '(TM10X)', sep = ' '))
genecounts <- data.frame(table(tm10x.markers$gene))
nonreduns <- subset(genecounts, Freq == 1); head(nonreduns); nrow(nonreduns) # 2201
tm10x.markers <- subset(tm10x.markers, gene %in% nonreduns$Var1)
tm10x.markers <- droplevels(tm10x.markers); nrow(tm10x.markers) # 2201

tmss2.markers <- read.delim('../../../mouse/tabulamuris_ss2/degs/markers.mouse.ss2.celltype_refined.MAST.txt')
tmss2.markers <- subset(tmss2.markers, p_val_adj <= 0.05)
tmss2.markers$cluster <- factor(tmss2.markers$cluster, levels = unique(tmss2.markers$cluster))
tmss2.markers$cluster <- mapvalues(tmss2.markers$cluster, from = levels(tmss2.markers$cluster), to = paste(levels(tmss2.markers$cluster), '(TMSS2)', sep = ' '))
genecounts <- data.frame(table(tmss2.markers$gene))
nonreduns <- subset(genecounts, Freq == 1); head(nonreduns); nrow(nonreduns) # 5013
tmss2.markers <- subset(tmss2.markers, gene %in% nonreduns$Var1)
tmss2.markers <- droplevels(tmss2.markers); nrow(tmss2.markers) # 5013


mouseObj.markers <- rbind(mca.markers, tm10x.markers, tmss2.markers)
levels(mouseObj.markers$cluster)
head(mouseObj.markers); nrow(mouseObj.markers) # 7937

conserved <- read.delim('../../../2.genelists/mouse/diopt_mousetohuman.usegenes.txt')
conserved <- conserved[, c(1,8)]
rownames(conserved) <- conserved[, 1]
head(conserved, n=3); nrow(conserved) # 15840
mouseObj.markers <- subset(mouseObj.markers, gene %in% conserved$MouseGeneID)
head(mouseObj.markers); nrow(mouseObj.markers) # 7367
summary(mouseObj.markers$cluster)
# Progenitors   Erythroids    Monocytes   Macrophage    pDC   Basophil    T cells   NK cells    B cells   Plasma cells
# 92            250           70          35            22    22          77        21          22        22
# Progenitors   Granulocytopoietic cells    Neutrophils   Erythroids    Monocytes   DC    pDC   Basophil    T cells   B cells   Plasma cells
# 188           47                          121           287           241         59    236   113         95        488       177
# Progenitors   Granulocytopoietic cells    Neutrophils   Monocytes    DC   pDC   Basophil  T cells   NK cells    B cells   Plasma cells
# 2385          322                         159           251          71   161   72        38        47          1122      54 
#write.table(mouseObj.markers, 'mouse_markers_celltype.orthologous.txt', row.names = F, col.names = T, sep = '\t', quote = F)


### Gene set ###
library(GSEABase)
library(GSVA)

gsInfo <- data.frame(row.names = levels(mouseObj.markers$cluster), Source = levels(mouseObj.markers$cluster)); head(gsInfo)
gsInfo <- AnnotatedDataFrame(data=gsInfo); gsInfo

targetgene_w <- data.frame(matrix(nrow = max(summary(mouseObj.markers$cluster)), ncol = length(levels(mouseObj.markers$cluster))))
colnames(targetgene_w) <- levels(mouseObj.markers$cluster)

gs_list <- list()
head(conserved, n=3)
for (i in c(1:length(levels(mouseObj.markers$cluster)))) {
  celltype <- levels(mouseObj.markers$cluster)[i]
  refgene <- as.character(subset(mouseObj.markers, cluster == celltype)$gene)
  #targetgene <- as.character(subset(conserved, MouseGeneID %in% refgene)$HumanSymbol)
  targetgene <- as.character(conserved[refgene, 'HumanSymbol'])
  if (length(targetgene) > max(summary(mouseObj.markers$cluster))) {
    targetgene <- targetgene[1:max(summary(mouseObj.markers$cluster))]
  }
  print (length(targetgene))
  gs_list[[i]] <- GeneSet(targetgene, setName = celltype)
  if (length(targetgene) < max(summary(mouseObj.markers$cluster))) {
    targetgene <- c(targetgene, rep('', time = max(summary(mouseObj.markers$cluster))-length(targetgene)))
  }
  targetgene_w[, celltype] <- targetgene
}
head(targetgene_w)
#write.table(targetgene_w, 'mouse_markers_celltype.orthologous.HumanSymbol.txt', row.names = F, col.names = T, sep = '\t', quote = F)
mySignatures <- GeneSetCollection(object = gs_list) #


### GSVA ###
# HCA
humanexpr <- readRDS('../../../human/hca/tmp/expr.pseudobulk.Rds')
humanexpr <- humanexpr[intersect(as.character(conserved$HumanSymbol), rownames(humanexpr)), ]
head(humanexpr, n = 3); dim(humanexpr) # 15837    15
humanexpr$sums <- rowSums(humanexpr)
humanexpr <- subset(humanexpr, sums > 0 )
humanexpr$sums <- NULL
head(humanexpr, n = 3); dim(humanexpr) # 15812    15

label <- data.frame(row.names = colnames(humanexpr), celltype = colnames(humanexpr)); label
label <- AnnotatedDataFrame(data=label)

exprSet <- ExpressionSet(as.matrix(humanexpr), phenoData = label, annotation = 'Symbol'); exprSet

es.gsva <- gsva(expr = exprSet, gset.idx.list = mySignatures,
                method='gsva', verbose=T, kcdf = 'Gaussian', parallel.sz =1)
gsvaScores <- as.data.frame(t(exprs(es.gsva)))
head(gsvaScores)
gsvaScores_w <- data.frame(celltype = rownames(gsvaScores), gsvaScores, check.rows = F, check.names = F)
#write.table(gsvaScores_w, 'gsva_mousetohuman_celltype.hca.txt', quote = F, sep = '\t', row.names = F, col.names = T)
min(gsvaScores); max(gsvaScores)

pheatmap(gsvaScores, cluster_rows = F, cluster_cols = F, scale = 'none', treeheight_row = 10, treeheight_col = 10, gaps_col = c(11,22),
         cellwidth = 10, cellheight = 10, fontsize_row = 9, fontsize_col = 9, border_color = NA,
         clustering_distance_rows = 'euclidean', clustering_distance_cols = 'euclidean', clustering_method = 'complete', clustering_callback = callback, 
         breaks = c(-0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1, 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8),
         color = colorRampPalette(c("#3a5fcd", 'white', "#ee0000"))(16), filename = 'gsva_mousetohuman_celltype.hca.pdf')
dev.off()


# HCL
humanexpr <- readRDS('../../../human/hcl/tmp/celltype_refined.pseudo.Rds')
humanexpr <- humanexpr[intersect(as.character(conserved$HumanSymbol), rownames(humanexpr)), ]
head(humanexpr, n = 3); dim(humanexpr) # 15157    11
humanexpr$sums <- rowSums(humanexpr)
humanexpr <- subset(humanexpr, sums > 0 )
humanexpr$sums <- NULL
head(humanexpr, n = 3); dim(humanexpr) # 12673    11

label <- data.frame(row.names = colnames(humanexpr), celltype = colnames(humanexpr)); label
label <- AnnotatedDataFrame(data=label)

exprSet <- ExpressionSet(as.matrix(humanexpr), phenoData = label, annotation = 'Symbol'); exprSet

es.gsva <- gsva(expr = exprSet, gset.idx.list = mySignatures,
                method='gsva', verbose=T, kcdf = 'Gaussian', parallel.sz =1)
gsvaScores <- as.data.frame(t(exprs(es.gsva)))
head(gsvaScores)
gsvaScores_w <- data.frame(celltype = rownames(gsvaScores), gsvaScores, check.rows = F, check.names = F)
#write.table(gsvaScores_w, 'gsva_mousetohuman_celltype.hcl.txt', quote = F, sep = '\t', row.names = F, col.names = T)
min(gsvaScores); max(gsvaScores)

pheatmap(gsvaScores, cluster_rows = F, cluster_cols = F, scale = 'none', treeheight_row = 10, treeheight_col = 10, gaps_col = c(11,22),
         cellwidth = 10, cellheight = 10, fontsize_row = 9, fontsize_col = 9, border_color = NA,
         clustering_distance_rows = 'euclidean', clustering_distance_cols = 'euclidean', clustering_method = 'complete', clustering_callback = callback, 
         breaks = c(-0.9, -0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1, 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9),
         color = colorRampPalette(c("#3a5fcd", 'white', "#ee0000"))(18), filename = 'gsva_mousetohuman_celltype.hcl.pdf')
dev.off()


