library(Seurat)
library(ggplot2)
library(plyr)
library(dplyr)

mca <- readRDS('../mca/tmp/mca.harmony.clustered.v2.Rds')
#mca <- subset(mca, celltype_refined != 'Neutrophils')
mca@meta.data <- droplevels(mca@meta.data)
head(mca@meta.data)
DimPlot(mca, label = T)

tm_10x <- readRDS('../tabulamuris_10x/tmp/tm_10x.harmony.clustered.Rds')
head(tm_10x@meta.data)
DimPlot(tm_10x, label = T)

tm_ss2 <- readRDS('../tabulamuris_ss2/tmp/tm_ss2.harmony.clustered.Rds')
head(tm_ss2@meta.data)
DimPlot(tm_ss2, label = T)


### pseudo-bulk ###
# mca
mca_pseudobulk <- data.frame(matrix(nrow = nrow(mca), ncol = length(levels(Idents(mca)))))
colnames(mca_pseudobulk) <- levels(Idents(mca))
rownames(mca_pseudobulk) <- rownames(mca)
mca_pseudobulk[is.na(mca_pseudobulk)] <- 0
head(mca_pseudobulk)

exprs_mca <- data.frame(as.matrix(GetAssayData(mca, slot = 'data')), check.rows = F, check.names = F)
exprs_mca <- exp(exprs_mca)-1; head(colSums(exprs_mca))
exprs_mca[1:4, 1:4]

for (ct in levels(Idents(mca)) ) {
  bcs <- rownames( subset(mca@meta.data, celltype_refined == ct) )
  exprs_mca_tmp <- exprs_mca[, bcs]
  exprs_mca_tmp[is.na(exprs_mca_tmp)] <- 0
  mca_pseudobulk[, ct] <- rowMeans(exprs_mca_tmp)
  remove(exprs_mca_tmp)
}
colnames(mca_pseudobulk) <- paste(colnames(mca_pseudobulk), ' (MCA)', sep = '')
head(mca_pseudobulk); dim(mca_pseudobulk)
colSums(mca_pseudobulk)
saveRDS(mca_pseudobulk, '../mca/tmp/mca_pseudobulk.Rds')
#saveRDS(mca_pseudobulk, '../mca/tmp/mca_pseudobulk_noNeu.Rds')

# 10x
tm_10x_pseudobulk <- data.frame(matrix(nrow = nrow(tm_10x), ncol = length(levels(Idents(tm_10x)))))
colnames(tm_10x_pseudobulk) <- levels(Idents(tm_10x))
rownames(tm_10x_pseudobulk) <- rownames(tm_10x)
tm_10x_pseudobulk[is.na(tm_10x_pseudobulk)] <- 0
head(tm_10x_pseudobulk)

exprs_10x <- data.frame(as.matrix(GetAssayData(tm_10x, slot = 'data')), check.rows = F, check.names = F)
exprs_10x <- exp(exprs_10x)-1
exprs_10x[1:4, 1:4]

for (ct in levels(Idents(tm_10x)) ) {
  bcs <- rownames( subset(tm_10x@meta.data, celltype_refined == ct) )
  exprs_10x_tmp <- exprs_10x[, bcs]
  exprs_10x_tmp[is.na(exprs_10x_tmp)] <- 0
  tm_10x_pseudobulk[, ct] <- rowMeans(exprs_10x_tmp)
  remove(exprs_10x_tmp)
}
colnames(tm_10x_pseudobulk) <- paste(colnames(tm_10x_pseudobulk), ' (10X)', sep = '')
head(tm_10x_pseudobulk); dim(tm_10x_pseudobulk)
colSums(tm_10x_pseudobulk)
saveRDS(tm_10x_pseudobulk, '../tabulamuris_10x/tmp/tm_10x_pseudobulk.Rds')

# ss2
tm_ss2_pseudobulk <- data.frame(matrix(nrow = nrow(tm_ss2), ncol = length(levels(Idents(tm_ss2)))))
colnames(tm_ss2_pseudobulk) <- levels(Idents(tm_ss2))
rownames(tm_ss2_pseudobulk) <- rownames(tm_ss2)
tm_ss2_pseudobulk[is.na(tm_ss2_pseudobulk)] <- 0
head(tm_ss2_pseudobulk)

exprs_ss2 <- data.frame(as.matrix(GetAssayData(tm_ss2, slot = 'data')), check.rows = F, check.names = F)
exprs_ss2 <- exp(exprs_ss2)-1
exprs_ss2[1:4, 1:4]

for (ct in levels(Idents(tm_ss2)) ) {
  bcs <- rownames( subset(tm_ss2@meta.data, celltype_refined == ct) )
  exprs_ss2_tmp <- exprs_ss2[, bcs]
  exprs_ss2_tmp[is.na(exprs_ss2_tmp)] <- 0
  tm_ss2_pseudobulk[, ct] <- rowMeans(exprs_ss2_tmp)
  remove(exprs_ss2_tmp)
}
colnames(tm_ss2_pseudobulk) <- paste(colnames(tm_ss2_pseudobulk), ' (SS2)', sep = '')
head(tm_ss2_pseudobulk); dim(tm_ss2_pseudobulk)
colSums(tm_ss2_pseudobulk)
saveRDS(tm_ss2_pseudobulk, '../tabulamuris_ss2/tmp/tm_ss2_pseudobulk.Rds')

