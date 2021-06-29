library(Seurat)
library(ggplot2)
library(plyr)
library(dplyr)

dir.create('tmp')
dir.create('res')
dir.create('degs')

### hca ###
hca <- readRDS('../../../Cross-species/human-HCA/tmp/integrated.hca.alignment.Rds')
head(hca@meta.data)
hca@meta.data$celltype <- NULL
hca@meta.data$integrated_snn_res.0.3 <- NULL

Idents(hca) <- 'integrated_snn_res.0.7'

markers <- FindAllMarkers(object = hca, min.pct = 0.25, logfc.threshold = 0.25, only.pos = T, test.use = 'MAST') # 
write.table(markers, 'degs/integrated.markers.hca.res_0.7.MAST.txt', sep = '\t', quote = F, col.names = T, row.names = F) 
top10 <- markers %>% group_by(cluster) %>% top_n(10, avg_logFC)
DoHeatmap(object = hca, features = top10$gene, angle = 90, size = 2, raster = T, draw.lines = F)
ggsave('degs/integrated.markers.hca.res_0.7.MAST.pdf', units = 'cm', width = 50, height = 40)

#  0  1  2  3  4 
#  5  6  7  8  9
# 10 11 12 13 14
# 15 16 17 18 19
# 20 21 22 23 24
# 25 26 27 28 29
# 30 31

celltypes <- c('Monocytes', 'T cells', 'T cells', 'B cells', 'T cells',
               'CD8 T cells', 'CD8 T cells', 'NK cells', 'NK cells', 'Monocytes',
               'B cells', 'Granulocyte progenitor', 'T cells', 'pre-B cells', 'Erythroids',
               'Progenitors', 'DC', 'NK cells', 'Monocytes', 'NK cells',
               'pDC', 'pre-PC', 'Erythroids', 'Plasma cells', 'CD8 T cells',
               'pro-B cells', 'B cells', 'pro-B cells', 'Doublets-proB_Ery', 'Stromal cells',
               'Platelet', 'Doublets-Mono_B')

hca@meta.data$celltype_refined <- mapvalues(hca@meta.data$integrated_snn_res.0.7, from = c(0:31), to = celltypes)
hca@meta.data$celltype_refined <- factor(hca@meta.data$celltype_refined, 
                                         levels = c("Progenitors", "Granulocyte progenitor", "Erythroids", "Monocytes", "DC", "pDC", 
                                                    "T cells", "CD8 T cells", "NK cells", "pre-PC", "pro-B cells", "pre-B cells", "B cells", "Plasma cells", 
                                                    "Platelet", "Stromal cells", "Doublets-Mono_B", "Doublets-proB_Ery"))
levels(hca@meta.data$celltype_refined)
Idents(hca) <- 'celltype_refined'

hca <- subset(hca, celltype_refined %in% c("Progenitors", "Granulocyte progenitor", "Erythroids", "Monocytes", "DC", "pDC", 
                                           "T cells", "CD8 T cells", "NK cells", "pre-PC", "pro-B cells", "pre-B cells", "B cells", "Plasma cells", "Platelet"))
hca@meta.data <- droplevels(hca@meta.data)
summary(hca@meta.data$celltype_refined)
saveRDS(hca, 'tmp/integrated.hca.alignment.filtered.Rds')


DimPlot(hca, label = T, pt.size = 0.5, reduction = 'umap', label.size = 2)
ggsave('res/umap.1_2.celltype_refined.png', units = 'cm', width = 16, height = 10)
DimPlot(hca, label = T, dims = c(1,3), pt.size = 0.5, reduction = 'umap', label.size = 2)
ggsave('res/umap.1_3.celltype_refined.png', units = 'cm', width = 16, height = 10)
DimPlot(hca, label = T, dims = c(2,3), pt.size = 0.5, reduction = 'umap', label.size = 2)
ggsave('res/umap.2_3.celltype_refined.png', units = 'cm', width = 16, height = 10)


markers <- FindAllMarkers(object = hca, min.pct = 0.25, logfc.threshold = 0.25, only.pos = T, test.use = 'MAST') # 
write.table(markers, 'degs/integrated.markers.hca.celltype_refined.MAST.txt', sep = '\t', quote = F, col.names = T, row.names = F) 
top10 <- markers %>% group_by(cluster) %>% top_n(10, avg_logFC)
DoHeatmap(object = hca, features = top10$gene, angle = 90, size = 2, raster = T, draw.lines = F)
ggsave('degs/integrated.markers.hca.celltype_refined.MAST.png', units = 'cm', width = 50, height = 40)



### Cell counts ###
library(reshape2)
options(scipen = 100)
dir.create('stats')

cellcounts <- data.frame(celltype = names(summary(hca@meta.data$celltype_refined)), count = summary(hca@meta.data$celltype_refined))
head(cellcounts)
cellcounts <- melt(cellcounts)
cellcounts$celltype <- factor(cellcounts$celltype, levels = levels(hca@meta.data$celltype_refined))
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
#ggsave('stats/cellcounts_by_celltype_refined.pdf', units = 'cm', width = 7, height = 7)



### Pseudo-cell (10 cells) transformation ###
pseudocell_total <- data.frame(matrix(ncol = 0, nrow = length(rownames(hca))))
rownames(pseudocell_total) <- rownames(hca)

for (ct in levels(hca@meta.data$celltype_refined)) {
  tmpObj <- subset(hca, celltype_refined == ct)
  
  if (nrow(tmpObj@meta.data) > 30000) {
    tmpObj_exprs <- data.frame(matrix(ncol = 0, nrow = length(rownames(hca))))
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
  
  ct_pseudocells <- data.frame(matrix(ncol = rounds, nrow = length(rownames(hca))))
  colnames(ct_pseudocells) <- unlist(lapply(c(1:rounds), function (x) paste0(c(as.character(ct), x), collapse = ' pc') ))
  rownames(ct_pseudocells) <- rownames(hca)
  
  for (i in c(1:rounds)) {
    use_bc <- sample(bcs, size = 10)
    #tmp_exprs <- t(FetchData(hca, vars = rownames(hca), cells = use_bc))
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
tptk <- data.frame(matrix(ncol = 0, nrow = length(rownames(hca))))
rownames(tptk) <- rownames(hca)
for (i in c(1: trunc(nrow(hca@meta.data)/10000) )) {
  if (i == trunc(nrow(hca@meta.data)/10000)) {
    start <- 10000*(i-1)+1
    end <- nrow(hca@meta.data)
  } else {
    start <- 10000*(i-1)+1
    end <- 10000*i
  }
  tmp_bc <- rownames(hca@meta.data)[start:end]
  tmpObj_subset <- subset(hca, cells = tmp_bc)
  tmpObj_data <- exp(as.matrix(GetAssayData(tmpObj_subset, slot = 'data', assay = 'RNA'))) - 1
  tptk <- cbind(tptk, tmpObj_data)
  remove(tmpObj_subset)
  remove(tmpObj_data)
}
dim(tptk)
head(colSums(tptk[, c(1:10)]))
saveRDS(tptk, 'tmp/expr.Rds')


pseudobulk <- data.frame(matrix(ncol = 0, nrow = length(rownames(hca))))
rownames(pseudobulk) <- rownames(tptk)
for (ct in levels(hca@meta.data$celltype_refined) ) {
  bcs <- rownames( subset(hca@meta.data, celltype_refined == ct) )
  pseudoB <- rowMeans(tptk[, bcs])
  pseudobulk <- cbind(pseudobulk, pseudoB)
}
colnames(pseudobulk) <- levels(hca@meta.data$celltype_refined)
head(pseudobulk, n = 3); dim(pseudobulk)
saveRDS(pseudobulk, 'tmp/expr.pseudobulk.Rds')



library(Seurat)
library(ggplot2)
library(plyr)
library(dplyr)
library(RColorBrewer)

hca <- readRDS('tmp/integrated.hca.alignment.filtered.Rds')
head(hca@meta.data)

mycol <- colorRampPalette(brewer.pal(9, "Spectral"))(17)
mycol[-c(4, 6)]

# hcl
#   1       2           3     4       5       6     7     8       9         10    11        12        13        14        15        16    17
#         gran.prog    Ery    Neu    Mono    Mac    DC    pDC    T cell           NK cell                                B cell    PC

# hca mycol[-c(3, 6)]
#   1       2           3     4       5       6     7     8       9         10    11        12        13        14        15        16    17
# Prog    gran.prog    Ery           Mono           DC    pDC    T cell    CD8    NK cell    pre-PC    pro-B    pre-B    B cell    PC    platelet


DimPlot(hca, pt.size = 0.5, reduction = 'umap') + 
  scale_color_manual(values = mycol[-c(4, 6)] ) + 
  theme_void() + 
  theme(legend.position = 'none')
ggsave('res/umap.1_2.celltype_refined.void.png', units = 'cm', width = 8, height = 8)
DimPlot(hca, pt.size = .5, reduction = 'umap', label = T, label.size = 2) + 
  scale_color_manual(values = mycol[-c(4, 6)] )
ggsave('res/umap.1_2.celltype_refined.png', units = 'cm', width = 16, height = 10)




