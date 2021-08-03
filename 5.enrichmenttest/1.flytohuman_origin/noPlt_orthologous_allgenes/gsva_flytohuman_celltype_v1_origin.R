callback = function(hc, mat){
  sv = svd(t(mat))$v[,1]
  dend = reorder(as.dendrogram(hc), wts = sv)
  as.hclust(dend)
}

library(pheatmap)
library(plyr)
library(dplyr)

flyObj_lg.markers <- readRDS('../../celltype_v1_lg.Rds')
flyObj_lg.markers <- subset(flyObj_lg.markers, p_val_adj <= 0.05)
head(flyObj_lg.markers); nrow(flyObj_lg.markers) # 2193

counts <- data.frame(table(flyObj_lg.markers$gene))
reduns <- subset(counts, Freq > 1); head(reduns); nrow(reduns) # 455
flyObj_lg.markers <- subset(flyObj_lg.markers, !gene %in% as.character(reduns$Var1))
head(flyObj_lg.markers); nrow(flyObj_lg.markers) # 1108
flyObj_lg.markers$cluster <- paste(flyObj_lg.markers$cluster, 'LG', sep = ' ')


flyObj_circ.markers <- readRDS('../../celltype_v1_circ.Rds')
flyObj_circ.markers <- subset(flyObj_circ.markers, p_val_adj <= 0.05)
head(flyObj_circ.markers); nrow(flyObj_circ.markers) # 2632

counts <- data.frame(table(flyObj_circ.markers$gene))
reduns <- subset(counts, Freq > 1); head(reduns); nrow(reduns) # 450
flyObj_circ.markers <- subset(flyObj_circ.markers, !gene %in% as.character(reduns$Var1))
head(flyObj_circ.markers); nrow(flyObj_circ.markers) # 1632
flyObj_circ.markers$cluster <- paste(flyObj_circ.markers$cluster, 'CI', sep = ' ')

flyObj.markers <- rbind(flyObj_lg.markers, flyObj_circ.markers)
flyObj.markers$cluster <- factor(flyObj.markers$cluster, 
                                 levels = c(unique(flyObj_lg.markers$cluster),
                                            unique(flyObj_circ.markers$cluster)) )
head(flyObj.markers); nrow(flyObj.markers) # 2740


conserved <- read.delim('../../../2.genelists/fly/diopt_flytohuman.usegenes.txt')
conserved <- conserved[, c(1,8)]
rownames(conserved) <- conserved[, 1]
head(conserved); nrow(conserved) # 6385
flyObj.markers <- subset(flyObj.markers, gene %in% conserved$FlyGeneID); nrow(flyObj.markers) # 1951
summary(flyObj.markers$cluster)
# PH 1 LG   PH LG   PM LG   PM 120 LG   LM LG   CC LG   GST-rich LG   Adipohemocyte LG 
# 238       29      27      68          157     72      73            86
# PH 1 CI   PH CI   PM CI   PM 120 CI   LM CI   CC CI   GST-rich CI 
# 412       119     55      58          351     117     89
#write.table(flyObj.markers, 'gsva_flytohuman_celltype_v1_origin.orthologous.txt', row.names = F, col.names = T, sep = '\t', quote = F)


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
  #targetgene <- as.character(subset(conserved, FlyGeneID %in% refgene)$HumanSymbol)
  targetgene <- as.character(conserved[refgene, 'HumanSymbol'])
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
#write.table(targetgene_w, 'gsva_flytohuman_celltype_v1_origin.orthologous.allgenes.HumanIDs.txt', row.names = F, col.names = T, sep = '\t', quote = F)
mySignatures <- GeneSetCollection(object = gs_list) #


### GSVA ###
humanexpr <- readRDS('../../../../Model_species/human/merged_commongenes/tmp/expr.pseudobulk.Rds')
humanexpr <- humanexpr[intersect(as.character(conserved$HumanSymbol), rownames(humanexpr)), ]
head(humanexpr, n = 3); dim(humanexpr) # 6385   17
humanexpr$Platelet <- NULL ###
humanexpr$sums <- rowSums(humanexpr)
humanexpr <- subset(humanexpr, sums > 0 )
humanexpr$sums <- NULL
head(humanexpr, n = 3); dim(humanexpr) # 6383   16

label <- data.frame(row.names = colnames(humanexpr), celltype = colnames(humanexpr)); label
label <- AnnotatedDataFrame(data=label)

exprSet <- ExpressionSet(as.matrix(humanexpr), phenoData = label, annotation = 'Symbol'); exprSet

es.gsva <- gsva(expr = exprSet, gset.idx.list = mySignatures,
                method='gsva', verbose=T, kcdf = 'Gaussian', parallel.sz =1)
gsvaScores <- as.data.frame(t(exprs(es.gsva)))
head(gsvaScores)
gsvaScores_w <- data.frame(celltype = rownames(gsvaScores), gsvaScores, check.rows = F, check.names = F)
#write.table(gsvaScores_w, 'gsva_flytohuman_celltype.txt', quote = F, sep = '\t', row.names = F, col.names = T)
min(gsvaScores); max(gsvaScores)

pheatmap(gsvaScores, cluster_rows = F, cluster_cols = F, scale = 'none', treeheight_row = 10, treeheight_col = 10,
         cellwidth = 10, cellheight = 10, fontsize_row = 9, fontsize_col = 9, border_color = NA,
         clustering_distance_rows = 'euclidean', clustering_distance_cols = 'euclidean', clustering_method = 'complete', clustering_callback = callback, 
         breaks = c(-0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1, 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7),
         color = colorRampPalette(c("#3a5fcd", 'white', "#ee0000"))(14), filename = 'gsva_flytohuman_celltype.pdf')
dev.off()


