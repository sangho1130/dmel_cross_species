library(Seurat)
library(ggplot2)
library(plyr)
library(dplyr)

hcl <- readRDS('../hcl/tmp/human_hcl.harmony.Rds')
head(hcl@meta.data)
DimPlot(hcl, label = T)

hcl_pseudobulk <- data.frame(matrix(nrow = nrow(hcl), ncol = length(levels(Idents(hcl)))))
colnames(hcl_pseudobulk) <- levels(Idents(hcl))
rownames(hcl_pseudobulk) <- rownames(hcl)
hcl_pseudobulk[is.na(hcl_pseudobulk)] <- 0
head(hcl_pseudobulk); dim(hcl_pseudobulk)

exprs_hcl <- data.frame(as.matrix(GetAssayData(hcl, slot = 'data')), check.rows = F, check.names = F)
exprs_hcl <- as.matrix(GetAssayData(hcl, slot = 'data'))
exprs_hcl <- exp(exprs_hcl)-1; head(colSums(exprs_hcl))
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

saveRDS(hcl_pseudobulk, '../hcl/tmp/celltype_refined.pseudo.Rds')
