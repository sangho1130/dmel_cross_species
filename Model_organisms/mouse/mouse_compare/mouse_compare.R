library(Seurat)
library(ggplot2)
library(plyr)
library(dplyr)

mca <- readRDS('../mca/tmp/mca.harmony.clustered.v2.Rds')
#mca <- subset(mca, celltype_refined != 'Neutrophils')
#mca@meta.data <- droplevels(mca@meta.data)
head(mca@meta.data)
DimPlot(mca, label = T)

tm_10x <- readRDS('../tabulamuris_10x/tmp/tm_10x.harmony.clustered.Rds')
head(tm_10x@meta.data)
DimPlot(tm_10x, label = T)

tm_ss2 <- readRDS('../tabulamuris_ss2/tmp/tm_ss2.harmony.clustered.Rds')
head(tm_ss2@meta.data)
DimPlot(tm_ss2, label = T)


### common variable genes ###
commongenes <- intersect(rownames(mca), intersect(rownames(tm_10x), rownames(tm_ss2)))
vst_table <- data.frame(matrix(nrow = length(commongenes), ncol = 4))
rownames(vst_table) <- commongenes
colnames(vst_table) <- c('mca', 'tenx', 'ss2', 'merged')

vst_table$mca <- mca@assays$RNA@meta.features[commongenes, 'vst.variance.standardized']
vst_table$tenx <- tm_10x@assays$RNA@meta.features[commongenes, 'vst.variance.standardized']
vst_table$ss2 <- tm_ss2@assays$RNA@meta.features[commongenes, 'vst.variance.standardized']
#vst_table$merged <- rowMeans(vst_table[, c(1:3)])
vst_table$merged <- rowSums(vst_table[, c(1:3)])
vst_table <- vst_table[order(vst_table$merged, decreasing = T),]
head(vst_table)
varigenes <- rownames(vst_table)[1:2000]


### pseudo-bulk ###
# mca
mca_pseudobulk <- data.frame(matrix(nrow = 2000, ncol = length(levels(Idents(mca)))))
colnames(mca_pseudobulk) <- levels(Idents(mca))
rownames(mca_pseudobulk) <- varigenes
mca_pseudobulk[is.na(mca_pseudobulk)] <- 0
head(mca_pseudobulk)

exprs_mca <- data.frame(as.matrix(GetAssayData(mca, slot = 'data')), check.rows = F, check.names = F)
exprs_mca <- exp(exprs_mca)-1; head(colSums(exprs_mca))
exprs_mca <- exprs_mca[varigenes, ]
exprs_mca[1:4, 1:4]

for (ct in levels(Idents(mca)) ) {
  bcs <- rownames( subset(mca@meta.data, celltype_refined == ct) )
  exprs_mca_tmp <- exprs_mca[, bcs]
  exprs_mca_tmp[is.na(exprs_mca_tmp)] <- 0
  mca_pseudobulk[, ct] <- rowMeans(exprs_mca_tmp)
  remove(exprs_mca_tmp)
}
colnames(mca_pseudobulk) <- paste(colnames(mca_pseudobulk), ' (MCA)', sep = '')
head(mca_pseudobulk)

# 10x
tm_10x_pseudobulk <- data.frame(matrix(nrow = 2000, ncol = length(levels(Idents(tm_10x)))))
colnames(tm_10x_pseudobulk) <- levels(Idents(tm_10x))
rownames(tm_10x_pseudobulk) <- varigenes
tm_10x_pseudobulk[is.na(tm_10x_pseudobulk)] <- 0
head(tm_10x_pseudobulk)

exprs_10x <- data.frame(as.matrix(GetAssayData(tm_10x, slot = 'data')), check.rows = F, check.names = F)
exprs_10x <- exp(exprs_10x)-1
exprs_10x <- exprs_10x[varigenes, ]
exprs_10x[1:4, 1:4]

for (ct in levels(Idents(tm_10x)) ) {
  bcs <- rownames( subset(tm_10x@meta.data, celltype_refined == ct) )
  exprs_10x_tmp <- exprs_10x[, bcs]
  exprs_10x_tmp[is.na(exprs_10x_tmp)] <- 0
  tm_10x_pseudobulk[, ct] <- rowMeans(exprs_10x_tmp)
  remove(exprs_10x_tmp)
}
colnames(tm_10x_pseudobulk) <- paste(colnames(tm_10x_pseudobulk), ' (10X)', sep = '')
head(tm_10x_pseudobulk)

# ss2
tm_ss2_pseudobulk <- data.frame(matrix(nrow = 2000, ncol = length(levels(Idents(tm_ss2)))))
colnames(tm_ss2_pseudobulk) <- levels(Idents(tm_ss2))
rownames(tm_ss2_pseudobulk) <- varigenes
tm_ss2_pseudobulk[is.na(tm_ss2_pseudobulk)] <- 0
head(tm_ss2_pseudobulk)

exprs_ss2 <- data.frame(as.matrix(GetAssayData(tm_ss2, slot = 'data')), check.rows = F, check.names = F)
exprs_ss2 <- exp(exprs_ss2)-1
exprs_ss2 <- exprs_ss2[varigenes, ]
exprs_ss2[1:4, 1:4]

for (ct in levels(Idents(tm_ss2)) ) {
  bcs <- rownames( subset(tm_ss2@meta.data, celltype_refined == ct) )
  exprs_ss2_tmp <- exprs_ss2[, bcs]
  exprs_ss2_tmp[is.na(exprs_ss2_tmp)] <- 0
  tm_ss2_pseudobulk[, ct] <- rowMeans(exprs_ss2_tmp)
  remove(exprs_ss2_tmp)
}
colnames(tm_ss2_pseudobulk) <- paste(colnames(tm_ss2_pseudobulk), ' (SS2)', sep = '')
head(tm_ss2_pseudobulk)


# corr 
library(pheatmap)

pseudo_corr <- cor(mca_pseudobulk, cbind(tm_10x_pseudobulk, tm_ss2_pseudobulk), method = 'spearman')
head(pseudo_corr)
max(pseudo_corr); min(pseudo_corr)
pheatmap(pseudo_corr, scale = 'none', gaps_col = 11,
         color = colorRampPalette(c("#3a5fcd", '#fefebd', "#ee0000"))(n = 12),
         breaks = c(0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0),
         cellwidth = 10, cellheight = 10, border_color = NA,
         cluster_rows = F, cluster_cols = F)#, filename = 'pseudobulk.spearman.mca_noNeu.pdf')


pseudo_corr <- cor(cbind(mca_pseudobulk, tm_10x_pseudobulk, tm_ss2_pseudobulk), cbind(mca_pseudobulk, tm_10x_pseudobulk, tm_ss2_pseudobulk), method = 'spearman')
pseudo_corr[upper.tri(pseudo_corr)] <- NA
head(pseudo_corr)

pheatmap(pseudo_corr, scale = 'none', gaps_col = c(11, 22), gaps_row = c(11, 22), 
         color = colorRampPalette(c("#3a5fcd", '#fefebd', "#ee0000"))(n = 12),
         breaks = c(0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0),
         cellwidth = 10, cellheight = 10, border_color = NA, na_col = "white",
         cluster_rows = F, cluster_cols = F)#, filename = 'pseudobulk.spearman.mca_noNeu.pdf')


