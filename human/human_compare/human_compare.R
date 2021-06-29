library(Seurat)
library(ggplot2)
library(plyr)
library(dplyr)

hcl <- readRDS('../hcl/tmp/human_hcl.harmony.Rds')
head(hcl@meta.data)
DimPlot(hcl, label = T)

hca_meta <- readRDS('../hca/tmp/rna.metafeatures.Rds')
head(hca_meta)


### common variable genes ###
commongenes <- intersect(rownames(hcl), rownames(hca_meta)); length(commongenes)
vst_table <- data.frame(matrix(nrow = length(commongenes), ncol = 3))
rownames(vst_table) <- commongenes
colnames(vst_table) <- c('hcl', 'hca', 'merged')

vst_table$hcl <- hcl@assays$RNA@meta.features[commongenes, 'vst.variance.standardized']
vst_table$hca <- hca_meta[commongenes, 'vst.variance.standardized']
#vst_table$merged <- rowMeans(vst_table[, c(1:3)])
vst_table$merged <- rowSums(vst_table[, c(1:2)])
vst_table <- vst_table[order(vst_table$merged, decreasing = T),]
head(vst_table)
varigenes <- rownames(vst_table)[1:2000]


### pseudo-bulk ###
# hcl
hcl_pseudobulk <- data.frame(matrix(nrow = 2000, ncol = length(levels(Idents(hcl)))))
colnames(hcl_pseudobulk) <- levels(Idents(hcl))
rownames(hcl_pseudobulk) <- varigenes
hcl_pseudobulk[is.na(hcl_pseudobulk)] <- 0
head(hcl_pseudobulk)

exprs_hcl <- data.frame(as.matrix(GetAssayData(hcl, slot = 'data')), check.rows = F, check.names = F)
exprs_hcl <- as.matrix(GetAssayData(hcl, slot = 'data'))
exprs_hcl <- exp(exprs_hcl)-1; head(colSums(exprs_hcl))
exprs_hcl <- exprs_hcl[varigenes, ]
exprs_hcl[1:4, 1:4]; dim(exprs_hcl)

for (ct in levels(Idents(hcl)) ) {
  bcs <- rownames( subset(hcl@meta.data, celltype_refined == ct) )
  exprs_hcl_tmp <- exprs_hcl[, bcs]
  exprs_hcl_tmp[is.na(exprs_hcl_tmp)] <- 0
  hcl_pseudobulk[, ct] <- rowMeans(exprs_hcl_tmp)
  remove(exprs_hcl_tmp)
}
colnames(hcl_pseudobulk) <- paste(colnames(hcl_pseudobulk), ' (hcl)', sep = '')
head(hcl_pseudobulk); dim(hcl_pseudobulk)

# hca
hca_pseudobulk <- readRDS('../hca/tmp/expr.pseudobulk.Rds')
colnames(hca_pseudobulk) <- paste(colnames(hca_pseudobulk), ' (hca)', sep = '')
hca_pseudobulk <- hca_pseudobulk[varigenes, ]
head(hca_pseudobulk); dim(hca_pseudobulk)


# corr 
library(pheatmap)
pseudo_corr <- cor(hcl_pseudobulk, hca_pseudobulk, method = 'spearman')
head(pseudo_corr)
min(pseudo_corr); max(pseudo_corr)

pheatmap(pseudo_corr, scale = 'none',
         color = colorRampPalette(c("#3a5fcd", '#fefebd', "#ee0000"))(n = 12),
         breaks = c(0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0),
         cellwidth = 10, cellheight = 10, border_color = NA,
         cluster_rows = F, cluster_cols = F)#, filename = 'pseudobulk.spearman.pdf')

