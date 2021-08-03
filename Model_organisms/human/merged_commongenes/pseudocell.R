library(Seurat)
library(ggplot2)
library(plyr)
library(dplyr)


### human ###
human <- readRDS('tmp/human.harmony_Library.flt.Rds')
human@meta.data <- droplevels(human@meta.data)
Idents(human) <- 'celltype_refined'
head(human@meta.data, n=3) # filtered = 262820

### Cell counts ###
library(reshape2)
options(scipen = 100)

cellcounts <- data.frame(celltype = names(summary(human@meta.data$celltype_refined)), count = summary(human@meta.data$celltype_refined))
head(cellcounts)
cellcounts <- melt(cellcounts)
cellcounts$celltype <- factor(cellcounts$celltype, levels = levels(human@meta.data$celltype_refined))
head(cellcounts)

ggplot(cellcounts, aes(celltype, value)) +
  geom_bar(stat = 'identity') +
  scale_y_continuous(trans='log10', breaks = c(0, 1, 10, 100, 1000, 10000, 100000)) +
  labs(x = '', y = 'Count') +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text = element_text(colour = 'black'),
        axis.ticks = element_line(colour = 'black'),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5),
        legend.position = 'none')
ggsave('stats/2-5.cellcounts_by_celltype_refined.pdf', units = 'cm', width = 7, height = 7)



### Pseudo-cell (10 cells) transformation ###
pseudocell_total <- data.frame(matrix(ncol = 0, nrow = length(rownames(human))))
rownames(pseudocell_total) <- rownames(human)

for (ct in levels(human@meta.data$celltype_refined)) {
  tmpObj <- subset(human, celltype_refined == ct)
  
  if (nrow(tmpObj@meta.data) > 30000) {
    tmpObj_exprs <- data.frame(matrix(ncol = 0, nrow = length(rownames(human))))
    for (i in c(1: trunc(nrow(tmpObj@meta.data)/10000) )) {
      if (i == trunc(nrow(tmpObj@meta.data)/10000)) {
        start <- 10000*(i-1)+1
        end <- nrow(tmpObj@meta.data)
      } else {
        start <- 10000*(i-1)+1
        end <- 10000*i
      }
      tmp_bc <- rownames(tmpObj@meta.data)[start:end]
      tmpObj_subset <- subset(tmpObj, cells = tmp_bc)
      tmpObj_exprs <- cbind(tmpObj_exprs, as.matrix(GetAssayData(tmpObj_subset, slot = 'data', assay = 'RNA')))
      remove(tmpObj_subset)
    }
  } else {
    tmpObj_exprs <- as.matrix(GetAssayData(tmpObj, slot = 'data', assay = 'RNA'))
  }
  
  bcs <- rownames(tmpObj@meta.data)
  rounds <- trunc(length(bcs)/10)
  
  ct_pseudocells <- data.frame(matrix(ncol = rounds, nrow = length(rownames(human))))
  colnames(ct_pseudocells) <- unlist(lapply(c(1:rounds), function (x) paste0(c(as.character(ct), x), collapse = ' pc') ))
  rownames(ct_pseudocells) <- rownames(human)
  
  for (i in c(1:rounds)) {
    use_bc <- sample(bcs, size = 10)
    #tmp_exprs <- t(FetchData(human, vars = rownames(human), cells = use_bc))
    tmp_exprs <- tmpObj_exprs[, use_bc]
    tmp_exprs <- exp(tmp_exprs)-1; colSums(tmp_exprs)
    tmp_exprs <- rowSums(tmp_exprs); sum(tmp_exprs)
    ct_pseudocells[, i] <- tmp_exprs
    
    bcs <- setdiff(bcs, use_bc)
  }
  
  pseudocell_total <- cbind(pseudocell_total, ct_pseudocells)
  
  remove(tmpObj)
  remove(tmpObj_exprs)
}
dim(pseudocell_total)
summary(colSums(pseudocell_total))
saveRDS(pseudocell_total, 'tmp/celltype_refined.pseudocell10.Rds')


### expression matrix - TP10K ###
tptk <- data.frame(matrix(ncol = 0, nrow = length(rownames(human))))
rownames(tptk) <- rownames(human)
for (i in c(1: trunc(nrow(human@meta.data)/10000) )) {
  if (i == trunc(nrow(human@meta.data)/10000)) {
    start <- 10000*(i-1)+1
    end <- nrow(human@meta.data)
  } else {
    start <- 10000*(i-1)+1
    end <- 10000*i
  }
  tmp_bc <- rownames(human@meta.data)[start:end]
  tmpObj_subset <- subset(human, cells = tmp_bc)
  tmpObj_data <- exp(as.matrix(GetAssayData(tmpObj_subset, slot = 'data', assay = 'RNA'))) - 1
  tptk <- cbind(tptk, tmpObj_data)
  remove(tmpObj_subset)
  remove(tmpObj_data)
}
dim(tptk)
head(colSums(tptk[, c(1:10)]))
saveRDS(tptk, 'tmp/expr.Rds')


pseudobulk <- data.frame(matrix(ncol = 0, nrow = length(rownames(human))))
rownames(pseudobulk) <- rownames(tptk)
for (ct in levels(human@meta.data$celltype_refined) ) {
  bcs <- rownames( subset(human@meta.data, celltype_refined == ct) )
  pseudoB <- rowMeans(tptk[, bcs])
  pseudobulk <- cbind(pseudobulk, pseudoB)
}
colnames(pseudobulk) <- levels(human@meta.data$celltype_refined)
head(pseudobulk, n = 3); dim(pseudobulk)
saveRDS(pseudobulk, 'tmp/expr.pseudobulk.Rds')


###
library(Seurat)
library(ggplot2)
library(plyr)
library(dplyr)

### human ###
human <- readRDS('tmp/human.harmony_Library.flt.Rds')
human@meta.data <- droplevels(human@meta.data)
Idents(human) <- 'celltype_refined'
head(human@meta.data, n=3) # filtered = 262820

### Pseudo-cell (10 cells) transformation ###
pseudocell_total <- data.frame(matrix(ncol = 0, nrow = length(rownames(human))))
rownames(pseudocell_total) <- rownames(human)

for (ct in levels(human@meta.data$celltype_refined)) {
  tmpObj <- subset(human, celltype_refined == ct)
  tmpObj_exprs <- as.matrix(GetAssayData(tmpObj, slot = 'data', assay = 'RNA'))
  
  bcs <- rownames(tmpObj@meta.data)
  rounds <- trunc(length(bcs)/10)
  
  ct_pseudocells <- data.frame(matrix(ncol = rounds, nrow = length(rownames(human))))
  colnames(ct_pseudocells) <- unlist(lapply(c(1:rounds), function (x) paste0(c(as.character(ct), x), collapse = ' pc') ))
  rownames(ct_pseudocells) <- rownames(human)
  
  for (i in c(1:rounds)) {
    use_bc <- sample(bcs, size = 10)
    tmp_exprs <- tmpObj_exprs[, use_bc]
    tmp_exprs <- exp(tmp_exprs)-1; colSums(tmp_exprs)
    tmp_exprs <- rowSums(tmp_exprs); sum(tmp_exprs)
    ct_pseudocells[, i] <- tmp_exprs
    
    bcs <- setdiff(bcs, use_bc)
  }
  
  pseudocell_total <- cbind(pseudocell_total, ct_pseudocells)
  
  remove(tmpObj)
  remove(tmpObj_exprs)
}
head(pseudocell_total, n=3); dim(pseudocell_total)
summary(colSums(pseudocell_total))
pseudocell_total$sums <- rowSums(pseudocell_total)
pseudocell_total <- subset(pseudocell_total, sums != 0)
pseudocell_total$sums <- NULL
dim(pseudocell_total)
saveRDS(pseudocell_total, 'tmp/celltype_refined.pseudocell10.Rds')


