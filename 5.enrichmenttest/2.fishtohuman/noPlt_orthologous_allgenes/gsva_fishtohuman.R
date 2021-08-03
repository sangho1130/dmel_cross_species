callback = function(hc, mat){
  sv = svd(t(mat))$v[,1]
  dend = reorder(as.dendrogram(hc), wts = sv)
  as.hclust(dend)
}

library(pheatmap)
library(plyr)
library(dplyr)

fishObj.markers <- read.delim('../../../../Model_species/zebrafish/Tang_Q_et_al/InDrops/degs/celltype.txt')
fishObj.markers <- subset(fishObj.markers, p_val_adj <= 0.05)
fishObj.markers$cluster <- factor(fishObj.markers$cluster, levels = c("HSCs", "Erythroid", "Neutrophil", "Macrophage", "NK/T cell", "B cell"))
head(fishObj.markers); nrow(fishObj.markers) # 1782

counts <- data.frame(table(fishObj.markers$gene))
nonreduns <- subset(counts, Freq == 1); head(nonreduns); nrow(nonreduns) # 1144
fishObj.markers <- subset(fishObj.markers, gene %in% as.character(nonreduns$Var1))
head(fishObj.markers); nrow(fishObj.markers) # 1144

conserved <- read.delim('../../../2.genelists/fish/diopt_fishtohuman.usegenes.txt')
conserved <- conserved[, c(1,8)]
rownames(conserved) <- conserved[, 1]
head(conserved); nrow(conserved) # 12044
fishObj.markers <- subset(fishObj.markers, gene %in% conserved$ZebrafishGeneID); nrow(fishObj.markers) # 776
head(fishObj.markers); summary(fishObj.markers$cluster)
# HSCs    Erythroid   Neutrophil  Macrophage    NK/T cell   B cell 
# 43        123        265        288         30         27 
#write.table(fishObj.markers, 'gsva_fishtohuman_celltype.orthologous.txt', row.names = F, col.names = T, sep = '\t', quote = F)

### Gene set ###
library(GSEABase)
library(GSVA)

gsInfo <- data.frame(row.names = levels(fishObj.markers$cluster), Source = levels(fishObj.markers$cluster)); head(gsInfo)
gsInfo <- AnnotatedDataFrame(data=gsInfo); gsInfo

targetgene_w <- data.frame(matrix(nrow = max(summary(fishObj.markers$cluster)), ncol = length(levels(fishObj.markers$cluster))))
head(targetgene_w)
colnames(targetgene_w) <- levels(fishObj.markers$cluster)

gs_list <- list()
head(conserved, n=3)
for (i in c(1:length(levels(fishObj.markers$cluster)))) {
  
  celltype <- levels(fishObj.markers$cluster)[i]
  refgene <- as.character(subset(fishObj.markers, cluster == celltype)$gene)
  #targetgene <- as.character(subset(conserved, ZebrafishGeneID %in% refgene)$HumanSymbol)
  targetgene <- as.character(conserved[refgene, 'HumanSymbol'])
  if (length(targetgene) > max(summary(fishObj.markers$cluster))) {
    targetgene <- targetgene[1:max(summary(fishObj.markers$cluster))]
  }
  print (length(targetgene))
  gs_list[[i]] <- GeneSet(targetgene, setName = celltype)
  if (length(targetgene) < max(summary(fishObj.markers$cluster))) {
    targetgene <- c(targetgene, rep('', times = max(summary(fishObj.markers$cluster))-length(targetgene)))
  }
  targetgene_w[, celltype] <- targetgene
}
head(targetgene_w)
#write.table(targetgene_w, 'gsva_fishtohuman_celltype.orthologous.HumanSymbol.txt', row.names = F, col.names = T, sep = '\t', quote = F)
mySignatures <- GeneSetCollection(object = gs_list) #


### GSVA ###
humanexpr <- readRDS('../../../../Cross-species_V2/human/merged_commongenes/tmp/expr.pseudobulk.Rds')
humanexpr <- humanexpr[intersect(as.character(conserved$HumanSymbol), rownames(humanexpr)), ]
head(humanexpr, n = 3); dim(humanexpr) # 12044    17
humanexpr$Platelet <- NULL ###
humanexpr$sums <- rowSums(humanexpr)
humanexpr <- subset(humanexpr, sums > 0 )
humanexpr$sums <- NULL
head(humanexpr, n = 3); dim(humanexpr) # 12043    16

label <- data.frame(row.names = colnames(humanexpr), celltype = colnames(humanexpr)); label
label <- AnnotatedDataFrame(data=label)

exprSet <- ExpressionSet(as.matrix(humanexpr), phenoData = label, annotation = 'Symbol'); exprSet

es.gsva <- gsva(expr = exprSet, gset.idx.list = mySignatures,
                method='gsva', verbose=T, kcdf = 'Gaussian', parallel.sz =1)
gsvaScores <- as.data.frame(t(exprs(es.gsva)))
head(gsvaScores)
gsvaScores_w <- data.frame(celltype = rownames(gsvaScores), gsvaScores, check.rows = F, check.names = F)
#write.table(gsvaScores_w, 'gsva_fishtohuman_celltype.txt', quote = F, sep = '\t', row.names = F, col.names = T)
min(gsvaScores); max(gsvaScores)

pheatmap(gsvaScores, cluster_rows = F, cluster_cols = F, scale = 'none', treeheight_row = 10, treeheight_col = 10,
         cellwidth = 10, cellheight = 10, fontsize_row = 9, fontsize_col = 9, border_color = NA,
         clustering_distance_rows = 'euclidean', clustering_distance_cols = 'euclidean', clustering_method = 'complete', clustering_callback = callback, 
         breaks = c(-0.5, -0.4, -0.3, -0.2, -0.1, 0.0, 0.1, 0.2, 0.3, 0.4, 0.5),
         color = colorRampPalette(c("#3a5fcd", 'white', "#ee0000"))(10), filename = 'gsva_fishtohuman_celltype.pdf')


