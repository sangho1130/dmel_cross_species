library(Seurat)
library(ggplot2)
library(plyr)
library(dplyr)

mouse_flt <- readRDS('tmp/mouse.harmony.flt.Rds')
head(mouse_flt@meta.data, n=3)
DimPlot(mouse_flt, label = T)

### pseudo-bulk ###
mouse_pseudobulk <- data.frame(matrix(nrow = nrow(mouse_flt), ncol = length(levels(Idents(mouse_flt)))))
colnames(mouse_pseudobulk) <- levels(Idents(mouse_flt))
rownames(mouse_pseudobulk) <- rownames(mouse_flt)
mouse_pseudobulk[is.na(mouse_pseudobulk)] <- 0
head(mouse_pseudobulk)


exprs_mouse <- data.frame(as.matrix(GetAssayData(mouse_flt, slot = 'data')), check.rows = F, check.names = F)
dim(exprs_mouse)
exprs_mouse <- exp(exprs_mouse)-1
exprs_mouse[1:4, 1:4]
colSums(exprs_mouse[, c(1:4)])

for (ct in levels(Idents(mouse_flt)) ) {
  
  bcs_dp <- rownames( subset(mouse_flt@meta.data, celltype_refined == ct & orig.ident != 'TMSS2') )
  exprs_mouse_dp_tmp <- exprs_mouse[, bcs_dp]
  #colSums(exprs_mouse_dp_tmp)
  
  bcs_ss <- rownames( subset(mouse_flt@meta.data, celltype_refined == ct & orig.ident == 'TMSS2') )
  if (length(bcs_ss) != 0) {
    exprs_mouse_ss_tmp <- exprs_mouse[, bcs_ss]
    #colSums(exprs_mouse_ss_tmp)
    
    mouse_pseudobulk[, ct] <- (rowMeans(exprs_mouse_dp_tmp) + rowMeans(exprs_mouse_ss_tmp))/2
    remove(exprs_mouse_dp_tmp, exprs_mouse_ss_tmp)
  } else {
    mouse_pseudobulk[, ct] <- rowMeans(exprs_mouse_dp_tmp)
    remove(exprs_mouse_dp_tmp)
  }
}
head(mouse_pseudobulk); dim(mouse_pseudobulk)
colSums(mouse_pseudobulk)
saveRDS(mouse_pseudobulk, 'tmp/mouse_pseudobulk_v2.Rds')


###
### Pseudo-cell (10 cells) transformation ###
pseudocell_total <- data.frame(matrix(ncol = 0, nrow = length(rownames(mouse_flt))))
rownames(pseudocell_total) <- rownames(mouse_flt)
pseudocell_total

for (ct in levels(mouse_flt@meta.data$celltype_refined)) {
  tmpObj <- subset(mouse_flt, celltype_refined == ct)
  tmpObj_exprs <- as.matrix(GetAssayData(tmpObj, slot = 'data', assay = 'RNA'))
  
  
  ### droplet
  bcs <- rownames(subset(tmpObj@meta.data, orig.ident != 'TMSS2'))
  rounds <- trunc(length(bcs)/10)
  
  ct_pseudocells <- data.frame(matrix(ncol = rounds, nrow = length(rownames(mouse_flt))))
  colnames(ct_pseudocells) <- unlist(lapply(c(1:rounds), function (x) paste0(c(as.character(ct), x), collapse = ' pc') ))
  rownames(ct_pseudocells) <- rownames(mouse_flt)
  
  lastnum <- 0
  for (i in c(1:rounds)) {
    use_bc <- sample(bcs, size = 10)
    tmp_exprs <- tmpObj_exprs[, use_bc]
    tmp_exprs <- exp(tmp_exprs)-1; colSums(tmp_exprs)
    tmp_exprs <- rowSums(tmp_exprs); sum(tmp_exprs)
    ct_pseudocells[, i] <- tmp_exprs
    
    bcs <- setdiff(bcs, use_bc)
    lastnum <- i
  }

  pseudocell_total <- cbind(pseudocell_total, ct_pseudocells)
  
  ### SS2
  bcs <- rownames(subset(tmpObj@meta.data, orig.ident == 'TMSS2'))
  rounds <- trunc(length(bcs)/10)
  if (length(bcs) != 0) {
    if (rounds > 0) {
      ct_pseudocells <- data.frame(matrix(ncol = rounds, nrow = length(rownames(mouse_flt))))
      colnames(ct_pseudocells) <- unlist(lapply(c(1:rounds), function (x) paste0(c(as.character(ct), x + lastnum), collapse = ' pc') ))
      rownames(ct_pseudocells) <- rownames(mouse_flt)
      
      for (i in c(1:rounds)) {
        use_bc <- sample(bcs, size = 10)
        tmp_exprs <- tmpObj_exprs[, use_bc]
        tmp_exprs <- exp(tmp_exprs)-1; colSums(tmp_exprs)
        tmp_exprs <- rowSums(tmp_exprs); sum(tmp_exprs)
        ct_pseudocells[, i] <- tmp_exprs
        
        bcs <- setdiff(bcs, use_bc)
      }
      
      pseudocell_total <- cbind(pseudocell_total, ct_pseudocells)
    }
  }
  
  remove(tmpObj)
  remove(tmpObj_exprs)
}
pseudocell_total[1:4, 1:10]; dim(pseudocell_total)
summary(colSums(pseudocell_total))
saveRDS(pseudocell_total, 'tmp/celltype_refined.pseudocell10_v2.Rds')


