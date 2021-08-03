library(Seurat)
library(plyr)
library(reshape2)
library(ggplot2)

dir.create('tmp')
dir.create('stats')

myObj <- readRDS('../../Drop-seq_merge_harmony_6.22/tmp/dropseq.combined.harmony_flt.Rds')
median(myObj@meta.data$nCount_RNA) # 5761 UMIs per cell
median(myObj@meta.data$nFeature_RNA) # 1476 genes per cell

head(myObj@meta.data)
Idents(myObj) <- 'celltype'
summary(Idents(myObj))

bloodcells <- subset(myObj, celltype != 'PSC'); bloodcells@meta.data <- droplevels(bloodcells@meta.data)
remove(myObj)
summary(Idents(bloodcells))
head(bloodcells@meta.data)

### Metadata preparation ###
bloodcells@meta.data$celltype_v2 <- mapvalues(bloodcells@meta.data$subclustering,
                                              from = levels(bloodcells@meta.data$subclustering),
                                              to = c("PH 1", "PH", "PH", "PH", "PH", "PH", "PM", "PM 120", "PM 120", "PM 120", 
                                                     "LM", "LM", "CC", "CC", "GST-rich", "Adipohemocyte"))

bloodcells@meta.data$celltype_v3 <- bloodcells@meta.data$celltype_v2
bloodcells@meta.data$origin_v2 <- mapvalues(bloodcells@meta.data$origin, from = levels(bloodcells@meta.data$origin), to = c("L", "C"))

bloodcells@meta.data$celltype_v3 <- unlist(lapply(c(1:nrow(bloodcells@meta.data)), 
                                                  function (x) paste0(c(as.character(bloodcells@meta.data$celltype_v2[x]), ' (', as.character(bloodcells@meta.data$origin_v2[x]), ')'), collapse = '') ))
bloodcells@meta.data$celltype_v3 <- factor(bloodcells@meta.data$celltype_v3,
                                           levels = c("PH 1 (L)", "PH (L)", "PM (L)", "PM 120 (L)", "LM (L)", "CC (L)", "GST-rich (L)", "Adipohemocyte (L)",
                                                      "PH 1 (C)", "PH (C)", "PM (C)", "PM 120 (C)", "LM (C)", "CC (C)", "GST-rich (C)"))
levels(bloodcells@meta.data$celltype_v3)

### Cell counts ###
cellcounts <- data.frame(celltype = names(summary(bloodcells@meta.data$celltype_v3)), count = summary(bloodcells@meta.data$celltype_v3))
cellcounts <- melt(cellcounts)
cellcounts$celltype <- factor(cellcounts$celltype, levels = levels(bloodcells@meta.data$celltype_v3))
cellcounts$variable <- mapvalues(cellcounts$celltype, from = levels(cellcounts$celltype), to = c('L', 'L', 'L', 'L', 'L', 'L', 'L', 'L', 'C', 'C', 'C', 'C', 'C', 'C', 'C'))
head(cellcounts)

#summary(bloodcells@meta.data$celltype_v3)
ggplot(cellcounts, aes(celltype, value, fill = variable)) +
  geom_bar(stat = 'identity') +
  scale_fill_manual(values = c('#ffa500', '#7ac5cd')) +
  scale_y_continuous(trans='log10', breaks = c(0, 1, 10, 100, 1000, 10000)) +
  labs(x = '', y = 'Count') +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text = element_text(colour = 'black'),
        axis.ticks = element_line(colour = 'black'),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5),
        legend.position = 'none')
#ggsave('stats/cellcounts_by_celltype_origin_v2.pdf', units = 'cm', width = 7, height = 6)


### Pseudo-cell (10 cells) transformation ###
#expr_total <- as.matrix(GetAssayData(bloodcells, slot = 'data', assay = 'RNA')); dim(expr_total); expr_total[1:4, 1:4]
#expr_total <- exp(expr_total) -1
#saveRDS(expr_total, 'tmp/exprs.Rds')
expr_total <- readRDS('tmp/exprs.Rds')
expr_total[1:4, 1:4]; dim(expr_total)

pseudocell_total <- data.frame(matrix(ncol = 0, nrow = length(rownames(bloodcells))))
rownames(pseudocell_total) <- rownames(bloodcells)

for (ct in levels(bloodcells@meta.data$celltype_v2)) {
  bcs <- rownames(subset(bloodcells@meta.data, celltype_v2 == ct))
  rounds <- trunc(length(bcs)/10)
  
  ct_pseudocells <- data.frame(matrix(ncol = rounds, nrow = length(rownames(bloodcells))))
  colnames(ct_pseudocells) <- unlist(lapply(c(1:rounds), function (x) paste0(c(as.character(ct), x), collapse = ' pc') ))
  rownames(ct_pseudocells) <- rownames(bloodcells)
  
  for (i in c(1:rounds)) {
    use_bc <- sample(bcs, size = 10)
    tmp_exprs <- expr_total[, use_bc]
    tmp_exprs <- rowSums(tmp_exprs); sum(tmp_exprs)
    ct_pseudocells[, i] <- tmp_exprs
    
    bcs <- setdiff(bcs, use_bc)
  }
  
  pseudocell_total <- cbind(pseudocell_total, ct_pseudocells)
  remove(ct_pseudocells)
}
head(pseudocell_total)
dim(pseudocell_total)
colnames(pseudocell_total)
#saveRDS(pseudocell_total, 'tmp/celltype_v2_origin_pseudocell10.Rds')
#saveRDS(pseudocell_total, 'tmp/celltype_v2_pseudocell10.Rds')


