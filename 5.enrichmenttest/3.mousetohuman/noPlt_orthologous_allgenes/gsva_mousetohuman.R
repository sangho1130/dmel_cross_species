callback = function(hc, mat){
  sv = svd(t(mat))$v[,1]
  dend = reorder(as.dendrogram(hc), wts = sv)
  as.hclust(dend)
}

library(pheatmap)
library(plyr)
library(dplyr)

mouseObj.markers <- read.delim('../../../../Model_species/mouse/merged_commongenes_v3/degs/markers.mouse.celltype_refined.DBflt.txt')
mouseObj.markers$cluster <- factor(mouseObj.markers$cluster, levels = unique(mouseObj.markers$cluster))
levels(mouseObj.markers$cluster)

genecounts <- data.frame(table(mouseObj.markers$gene))
nonreduns <- subset(genecounts, Freq == 1); head(nonreduns); nrow(nonreduns) # 2107
mouseObj.markers <- subset(mouseObj.markers, gene %in% nonreduns$Var1)
mouseObj.markers <- droplevels(mouseObj.markers); nrow(mouseObj.markers) # 2107
head(mouseObj.markers); nrow(mouseObj.markers) # 2107

conserved <- read.delim('../../../2.genelists/mouse/diopt_mousetohuman.usegenes.txt')
conserved <- conserved[, c(1,8)]
rownames(conserved) <- conserved[, 1]
head(conserved, n=3); nrow(conserved) # 10380
mouseObj.markers <- subset(mouseObj.markers, gene %in% conserved$MouseGeneID)
head(mouseObj.markers); nrow(mouseObj.markers) # 1995
summary(mouseObj.markers$cluster)
# Progenitors   Granulocytopoietic cells    Neutrophils   Erythroids    Monocytes   Macrophage    DC    pDC 
# 590           76                          152           210           61          78            88    189
# Basophil    T cells   NK cells    B cells   Plasma cells 
# 52          64        47          265       123
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
humanexpr <- readRDS('../../../../Cross-species_V2/human/merged_commongenes/tmp/expr.pseudobulk.Rds')
humanexpr <- humanexpr[intersect(as.character(conserved$HumanSymbol), rownames(humanexpr)), ]
head(humanexpr, n = 3); dim(humanexpr) # 10380    17
humanexpr$Platelet <- NULL ###
humanexpr$sums <- rowSums(humanexpr)
humanexpr <- subset(humanexpr, sums > 0 )
humanexpr$sums <- NULL
head(humanexpr, n = 3); dim(humanexpr) # 10380    16

label <- data.frame(row.names = colnames(humanexpr), celltype = colnames(humanexpr)); label
label <- AnnotatedDataFrame(data=label)

exprSet <- ExpressionSet(as.matrix(humanexpr), phenoData = label, annotation = 'Symbol'); exprSet

es.gsva <- gsva(expr = exprSet, gset.idx.list = mySignatures,
                method='gsva', verbose=T, kcdf = 'Gaussian', parallel.sz =1)
gsvaScores <- as.data.frame(t(exprs(es.gsva)))
head(gsvaScores)
gsvaScores_w <- data.frame(celltype = rownames(gsvaScores), gsvaScores, check.rows = F, check.names = F)
#write.table(gsvaScores_w, 'gsva_mousetohuman_celltype.txt', quote = F, sep = '\t', row.names = F, col.names = T)
min(gsvaScores); max(gsvaScores)

pheatmap(gsvaScores, cluster_rows = F, cluster_cols = F, scale = 'none', treeheight_row = 10, treeheight_col = 10, 
         cellwidth = 10, cellheight = 10, fontsize_row = 9, fontsize_col = 9, border_color = NA,
         clustering_distance_rows = 'euclidean', clustering_distance_cols = 'euclidean', clustering_method = 'complete', clustering_callback = callback, 
         breaks = c(-0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1, 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7),
         color = colorRampPalette(c("#3a5fcd", 'white', "#ee0000"))(14), filename = 'gsva_mousetohuman_celltype.pdf')


