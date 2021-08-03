library(Seurat)
library(ggplot2)
library(plyr)
library(dplyr)
library(harmony)


### mca ###
mca_counts <- readRDS('/home/sangho/2020_newbuild/projects/drosophila/projectX_review/Cross-species/mouse-MCA/MCA1.1_adata.counts.PB_BM.v3.Rds')
mca_label <- readRDS('/home/sangho/2020_newbuild/projects/drosophila/projectX_review/Cross-species/mouse-MCA/MCA1.1_cell_info.usecells.Rds')
mca_label <- subset(mca_label, tissue != 'AdultBoneMarrowcKit')

mca_counts <- mca_counts[, as.character(rownames(mca_label))]
mca_counts$sums <- rowSums(mca_counts)
mca_counts <- subset(mca_counts, sums != 0)
mca_counts$sums <- NULL

mca <- CreateSeuratObject(counts = mca_counts, project = "MCA")
mca <- AddMetaData(object = mca, metadata = label[as.character(rownames(mca_label)), ])
mca <- NormalizeData(mca, normalization.method = "LogNormalize", scale.factor = 10000)
mca@meta.data <- droplevels(mca@meta.data)
head(mca@meta.data); summary(mca@meta.data)

dir.create('tmp')
saveRDS(mca, 'tmp/mca.Rds')


### Standard workflow ###
dir.create('stats')
dir.create('tsne')
dir.create('tsne/compare')
dir.create('umap')
dir.create('umap/compare')
dir.create('degs')

mca <- readRDS('tmp/mca.Rds')
head(mca@meta.data); nrow(mca@meta.data) # 16144 cells

mca <- FindVariableFeatures(object = mca, selection.method = "vst", nfeatures = 2000)
mca <- ScaleData(object = mca, vars.to.regress = c('nCount_RNA'))
mca <- RunPCA(object = mca, npcs = 60)
mca <- JackStraw(object = mca, num.replicate = 100, dims = 60)
mca <- ScoreJackStraw(object = mca, dims = 1:60)
JackStrawPlot(mca, dims = 1:60) # 35 PCs 
ggsave('stats/2-4.pca.JackStrawPlot.mitoCut.pdf', units = 'cm', width = 25, height = 12)
saveRDS(mca, 'tmp/mca_2.Rds')

mca <- RunHarmony(mca, group.by.vars = "tissue")
saveRDS(mca, 'tmp/mca.harmony.Rds')

#mca <- readRDS('tmp/mca.harmony.Rds')
### UMAP ###
for (seed in c(210516150:210516169)){
  mca <- RunUMAP(mca, dims=1:35, reduction = "harmony", reduction.key='UMAP', n.components=3, min.dist=0.3, seed.use = seed)
  DimPlot(mca, group.by = "celltype_merged", reduction = "umap", label = T, label.size = 2)
  ggsave(paste0(c('umap/compare/umap.1_2.markers.seed', seed, '.png'), collapse = ''), units = 'cm', width = 14, height = 10)
} # 210516152

### tSNE ###
for (seed in c(210516100:210516109)){
  mca <- RunTSNE(mca, dims=1:35, reduction = "harmony", reduction.key='tSNE', dim.embed=3, seed.use = seed)
  DimPlot(mca, group.by = "celltype_merged", reduction = "tsne", label = T, label.size = 2)
  ggsave(paste0(c('tsne/compare/tsne.1_2.markers.seed', seed, '.png'), collapse = ''), units = 'cm', width = 14, height = 10)
} # 210516102


mca <- RunUMAP(mca, dims=1:35, reduction = "harmony", reduction.key='UMAP', n.components=3, min.dist=0.3, seed.use = 210516152)
DimPlot(object=mca, group.by = "dataset", reduction = 'umap')
ggsave('umap/umap.1_2.dataset.seed210516152.pdf', units = 'cm', width = 13.5, height = 10)
DimPlot(object=mca, group.by = "tissue", reduction = 'umap')
ggsave('umap/umap.1_2.tissue.seed210516152.pdf', units = 'cm', width = 14.5, height = 10)
DimPlot(object=mca, group.by = "celltype_merged", reduction = 'umap')
ggsave('umap/umap.1_2.celltype_merged.seed210516152.pdf', units = 'cm', width = 14.5, height = 10)

mca <- RunTSNE(mca, dims=1:35, reduction = "harmony", reduction.key='tSNE', dim.embed=3, seed.use = 210516102)
DimPlot(object=mca, group.by = "dataset", reduction = 'tsne')
ggsave('tsne/tsne.1_2.dataset.seed210516102.pdf', units = 'cm', width = 13.5, height = 10)
DimPlot(object=mca, group.by = "tissue", reduction = 'tsne')
ggsave('tsne/tsne.1_2.tissue.seed210516102.pdf', units = 'cm', width = 14.5, height = 10)
DimPlot(object=mca, group.by = "celltype_merged", reduction = 'tsne')
ggsave('tsne/tsne.1_2.celltype_merged.seed210516102.pdf', units = 'cm', width = 14.5, height = 10)
saveRDS(mca, 'tmp/mca.harmony.Rds')


### clustering ###
dir.create('res')
dir.create('res/compare')
mca <- FindNeighbors(object = mca, dims = 1:35)
for (res in c(1:12)){
  res <- res/10
  mca <- FindClusters(object = mca, resolution = res)
  DimPlot(mca, label = T, pt.size = 0.5, reduction = 'umap')
  ggsave(paste(c('res/compare/umap.seed210516152.1_2.res_',res,'.pdf'), collapse = ''), units = 'cm', width = 13.5, height = 10)
  DimPlot(mca, label = T, dims = c(1,3), pt.size = 0.5, reduction = 'umap')
  ggsave(paste(c('res/compare/umap.seed210516152.1_3.res_',res,'.pdf'), collapse = ''), units = 'cm', width = 13.5, height = 10)
  DimPlot(mca, label = T, dims = c(2,3), pt.size = 0.5, reduction = 'umap')
  ggsave(paste(c('res/compare/umap.seed210516152.2_3.res_',res,'.pdf'), collapse = ''), units = 'cm', width = 13.5, height = 10)
  
  DimPlot(mca, label = T, pt.size = 0.5, reduction = 'tsne')
  ggsave(paste(c('res/compare/tsne.seed210516102.1_2.res_',res,'.pdf'), collapse = ''), units = 'cm', width = 13.5, height = 10)
  DimPlot(mca, label = T, dims = c(1,3), pt.size = 0.5, reduction = 'tsne')
  ggsave(paste(c('res/compare/tsne.seed210516102.1_3.res_',res,'.pdf'), collapse = ''), units = 'cm', width = 13.5, height = 10)
  DimPlot(mca, label = T, dims = c(2,3), pt.size = 0.5, reduction = 'tsne')
  ggsave(paste(c('res/compare/tsne.seed210516102.2_3.res_',res,'.pdf'), collapse = ''), units = 'cm', width = 13.5, height = 10)
}

mca <- readRDS('tmp/mca.harmony.Rds')
mca <- FindNeighbors(object = mca, dims = 1:35)

mca <- FindClusters(object = mca, resolution = 0.3)
DimPlot(mca, label = T, pt.size = 0.5, reduction = 'umap')
ggsave('res/umap.seed210516152.1_2.res_0.3.pdf', units = 'cm', width = 12.5, height = 10)
DimPlot(mca, label = T, dims = c(1,3), pt.size = 0.5, reduction = 'umap')
ggsave('res/umap.seed210516152.1_3.res_0.3.pdf', units = 'cm', width = 12.5, height = 10)
DimPlot(mca, label = T, dims = c(2,3), pt.size = 0.5, reduction = 'umap')
ggsave('res/umap.seed210516152.2_3.res_0.3.pdf', units = 'cm', width = 12.5, height = 10)

DimPlot(mca, label = T, pt.size = 0.5, reduction = 'tsne')
ggsave('res/tsne.seed210516102.1_2.res_0.3.pdf', units = 'cm', width = 12.5, height = 10)
DimPlot(mca, label = T, dims = c(1,3), pt.size = 0.5, reduction = 'tsne')
ggsave('res/tsne.seed210516102.1_3.res_0.3.pdf', units = 'cm', width = 12.5, height = 10)
DimPlot(mca, label = T, dims = c(2,3), pt.size = 0.5, reduction = 'tsne')
ggsave('res/tsne.seed210516102.2_3.res_0.3.pdf', units = 'cm', width = 12.5, height = 10)
saveRDS(mca, 'tmp/mca.harmony.Rds')


### markers ###
markers <- FindAllMarkers(object = mca, min.pct = 0.25, logfc.threshold = 0.25, only.pos = T)
markers %>% group_by(cluster) %>% top_n(2, avg_logFC)
write.table(markers, 'degs/markers.mca.res_0.3.txt', sep = '\t', quote = F, col.names = T, row.names = F) 
top10 <- markers %>% group_by(cluster) %>% top_n(10, avg_logFC)
DoHeatmap(object = mca, features = top10$gene, angle = 90, size = 5, raster = T, draw.lines = F)
ggsave('degs/markers.mca.res_0.3.pdf', units = 'cm', width = 30, height = 20)



library(Seurat)
library(ggplot2)
library(plyr)
library(dplyr)

mca <- readRDS('tmp/mca.harmony.Rds')
mca <- readRDS('tmp/mca.harmony.clustered.Rds')

VlnPlot(mca, pt.size = 0, features = c('S100a9', 'S100a8', 'Mmp9', 'Mmp8')) # Neutrophils
VlnPlot(mca, pt.size = 0, features = c('Gzma', 'Klrb1c', 'Ncr1')) # NK cells
VlnPlot(mca, pt.size = 0, features = c('Prss34', 'Mcpt8', 'Ccl3', 'Cpa3', 'Ccl4', 'Fcer1a')) # Mast cells
VlnPlot(mca, pt.size = 0, features = c('Siglech', 'Ccr9', 'Bst2', 'Pacsin1', 'Tcf4')) # pDC
VlnPlot(mca, pt.size = 0, features = c('Csf1r', 'C1qa', 'Apoe', 'Fn1')) # Mono macrophage
VlnPlot(mca, pt.size = 0, features = c('Alas1','Anxa3', 'Elane', 'Gatm', 'Gm11505', 'Hk3', 'Hp', 'Igsf6')) # GMP
VlnPlot(mca, pt.size = 0, features = c('Cd79a', 'Cd79b', 'Mzb1', 'Jchain', 'Iglv1', 'Iglc1', 'Iglc2', 'Derl3', 'Tnfrsf17')) # GMP
VlnPlot(mca, pt.size = 0, features = c('Cd79a', 'Cd79b', 'Cd3d', 'Cd3e', 'Cd3g'), group.by = 'RNA_snn_res.0.3')
VlnPlot(mca, pt.size = 0, features = c('Fcnb', 'Orm1', 'Mogat2', 'Cd34', 'Igfbp4', 'Eltd1'), group.by = 'RNA_snn_res.0.3')
VlnPlot(mca, pt.size = 0, features = c('Csf1r', 'C1qa', 'Apoe', 'Fn1'), group.by = 'RNA_snn_res.0.3')


mca <- subset(mca, RNA_snn_res.0.3 != '10')
mca@meta.data <- droplevels(mca@meta.data)
### rename ### 
# 0 1 2 3 4 5
# 6 7 8 9 11
# 12 13 14 15 16 17
clusters <- c('Neutrophils', 'Neutrophils', 'T cells', 'Progenitors', 'Monocytes', 'Neutrophils',
              'B cells', 'Progenitors', 'B cells', 'Erythroids', 'Neutrophils',
              'NK cells', 'pDC', 'Neutrophils', 'Basophil', 'Plasma cells', 'Macrophage')
mca@meta.data$celltype_refined <- mapvalues(mca@meta.data$RNA_snn_res.0.3, from = levels(mca@meta.data$RNA_snn_res.0.3), to = clusters)
mca@meta.data$celltype_refined <- factor(mca@meta.data$celltype_refined, 
                                         levels = c("Progenitors", "Neutrophils", "Erythroids", "Monocytes", "Macrophage", "pDC", "Basophil", "T cells", "NK cells", "B cells", "Plasma cells"))
Idents(mca) <- 'celltype_refined'


DimPlot(mca, label = T, pt.size = 0.5, reduction = 'umap', label.size = 2)
ggsave('res/umap.seed210516152.1_2.celltype_refined.pdf', units = 'cm', width = 14, height = 10)
DimPlot(mca, label = T, dims = c(1,3), pt.size = 0.5, reduction = 'umap', label.size = 2)
ggsave('res/umap.seed210516152.1_3.celltype_refined.pdf', units = 'cm', width = 14, height = 10)
DimPlot(mca, label = T, dims = c(2,3), pt.size = 0.5, reduction = 'umap', label.size = 2)
ggsave('res/umap.seed210516152.2_3.celltype_refined.pdf', units = 'cm', width = 14, height = 10)

DimPlot(mca, label = T, pt.size = 0.5, reduction = 'tsne', label.size = 2)
ggsave('res/tsne.seed210516102.1_2.celltype_refined.pdf', units = 'cm', width = 14, height = 10)
DimPlot(mca, label = T, dims = c(1,3), pt.size = 0.5, reduction = 'tsne', label.size = 2)
ggsave('res/tsne.seed210516102.1_3.celltype_refined.pdf', units = 'cm', width = 14, height = 10)
DimPlot(mca, label = T, dims = c(2,3), pt.size = 0.5, reduction = 'tsne', label.size = 2)
ggsave('res/tsne.seed210516102.2_3.celltype_refined.pdf', units = 'cm', width = 14, height = 10)
#saveRDS(mca, 'tmp/mca.harmony.clustered.v2.Rds')

library(Seurat)
library(plyr)
library(dplyr)
library(ggplot2)

mca <- readRDS('tmp/mca.harmony.clustered.v2.Rds')
levels(Idents(mca))
head(mca@meta.data)

markers <- FindAllMarkers(object = mca, min.pct = 0.25, logfc.threshold = 0.25, only.pos = T, test.use = 'MAST')
markers %>% group_by(cluster) %>% top_n(2, avg_logFC)
write.table(markers, 'degs/markers.mca.celltype_refined.MAST.txt', sep = '\t', quote = F, col.names = T, row.names = F) 
top10 <- markers %>% group_by(cluster) %>% top_n(10, avg_logFC)
DoHeatmap(object = mca, features = top10$gene, angle = 90, size = 3, raster = T, draw.lines = F)
ggsave('degs/markers.mca.celltype_refined.MAST.pdf', units = 'cm', width = 40, height = 30)



### Remove neutrophils ###
library(Seurat)
library(plyr)
library(dplyr)
library(ggplot2)

mca <- readRDS('tmp/mca.harmony.clustered.v2.Rds')
mca <- subset(mca, celltype_refined != 'Neutrophils')
mca@meta.data <- droplevels(mca@meta.data)
head(mca@meta.data)

DimPlot(mca, label = T, pt.size = 0.5, reduction = 'umap', label.size = 2)
ggsave('res/umap.seed210516152.1_2.celltype_refined_noNeu.pdf', units = 'cm', width = 13.5, height = 10)
DimPlot(mca, label = T, dims = c(1,3), pt.size = 0.5, reduction = 'umap', label.size = 2)
ggsave('res/umap.seed210516152.1_3.celltype_refined_noNeu.pdf', units = 'cm', width = 13.5, height = 10)
DimPlot(mca, label = T, dims = c(2,3), pt.size = 0.5, reduction = 'umap', label.size = 2)
ggsave('res/umap.seed210516152.2_3.celltype_refined_noNeu.pdf', units = 'cm', width = 13.5, height = 10)

DimPlot(mca, label = T, pt.size = 0.5, reduction = 'tsne', label.size = 2)
ggsave('res/tsne.seed210516102.1_2.celltype_refined_noNeu.pdf', units = 'cm', width = 13.5, height = 10)
DimPlot(mca, label = T, dims = c(1,3), pt.size = 0.5, reduction = 'tsne', label.size = 2)
ggsave('res/tsne.seed210516102.1_3.celltype_refined_noNeu.pdf', units = 'cm', width = 13.5, height = 10)
DimPlot(mca, label = T, dims = c(2,3), pt.size = 0.5, reduction = 'tsne', label.size = 2)
ggsave('res/tsne.seed210516102.2_3.celltype_refined_noNeu.pdf', units = 'cm', width = 13.5, height = 10)


levels(Idents(mca))
markers <- FindAllMarkers(object = mca, min.pct = 0.25, logfc.threshold = 0.25, only.pos = T, test.use = 'MAST')
markers %>% group_by(cluster) %>% top_n(2, avg_logFC)
write.table(markers, 'degs/markers.mca.celltype_refined_noNeu.MAST.txt', sep = '\t', quote = F, col.names = T, row.names = F) 
top10 <- markers %>% group_by(cluster) %>% top_n(10, avg_logFC)
DoHeatmap(object = mca, features = top10$gene, angle = 90, size = 3, raster = T, draw.lines = F)
ggsave('degs/markers.mca.celltype_refined_noNeu.MAST.pdf', units = 'cm', width = 35, height = 25)


### Pseudo-cell (10 cells) transformation ###
pseudocell_total <- data.frame(matrix(ncol = 0, nrow = length(rownames(mca))))
rownames(pseudocell_total) <- rownames(mca)

for (ct in levels(mca@meta.data$celltype_refined)) {
  tmpObj <- subset(mca, celltype_refined == ct)
  tmpObj_exprs <- as.matrix(GetAssayData(tmpObj, slot = 'data', assay = 'RNA'))
  
  bcs <- rownames(tmpObj@meta.data)
  rounds <- trunc(length(bcs)/10)
  
  ct_pseudocells <- data.frame(matrix(ncol = rounds, nrow = length(rownames(mca))))
  colnames(ct_pseudocells) <- unlist(lapply(c(1:rounds), function (x) paste0(c(as.character(ct), x), collapse = ' pc') ))
  rownames(ct_pseudocells) <- rownames(mca)
  
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
saveRDS(pseudocell_total, 'tmp/celltype_refined_noNeu.pseudocell10.Rds')




library(Seurat)
library(plyr)
library(dplyr)
library(ggplot2)
library(RColorBrewer)

mca <- readRDS('tmp/mca.harmony.clustered.v2.Rds')
mca <- subset(mca, celltype_refined != 'Neutrophils')
mca@meta.data <- droplevels(mca@meta.data)
head(mca@meta.data)

mca <- JackStraw(object = mca, num.replicate = 100, dims = 40)
mca <- ScoreJackStraw(object = mca, dims = 1:40)
JackStrawPlot(mca, dims = 1:40) # 23 PCs 
ggsave('stats/2-4.pca.JackStrawPlot.mitoCut.noNeu.pdf', units = 'cm', width = 20, height = 12)

mca <- RunUMAP(mca, dims=1:23, reduction = "harmony", reduction.key='UMAP', n.components=3, min.dist=0.3, seed.use = 210607109)
# 210607108, 210607109
saveRDS(mca, 'tmp/mca.harmony.clustered.v3.noNeu.Rds')


library(Seurat)
library(plyr)
library(dplyr)
library(ggplot2)
library(RColorBrewer)

mca <- readRDS('tmp/mca.harmony.clustered.v3.noNeu.Rds')
levels(Idents(mca))
#"Progenitors",                             "Erythroids",                "Monocytes", "Macrophage",       "pDC", "Basophil", "T cells", "NK cells", "B cells", "Plasma cells"
tm_merged <- readRDS('../tabulamuris_merged/tmp/tm_merged.harmony.Rds')
Idents(tm_merged) <- 'celltype_refined'
levels(Idents(tm_merged))
#"Progenitors", "Granulocytopoietic cells", "Erythroids", "Neutrophils", "Monocytes",               "DC", "pDC", "Basophil", "T cells", "NK cells", "B cells", "Plasma cells"

cols <- colorRampPalette(brewer.pal(9, "Spectral"))(13)
cols[-c(2,4,7)]

DimPlot(mca, pt.size = 0.5, reduction = 'umap', label = T, label.size = 2) + 
  scale_color_manual(values = cols[-c(2,4,7)])
ggsave('res/new.umap.seed210607109.1_2.celltype_refined_noNeu.pdf', units = 'cm', width = 13.5, height = 10)
DimPlot(mca, pt.size = 0.5, reduction = 'umap') + 
  scale_color_manual(values = cols[-c(2,4,7)]) + 
  theme_void() + 
  theme(legend.position = 'none')
ggsave('res/new.umap.seed210607109.1_2.celltype_refined_noNeu.void.png', units = 'cm', width = 8, height = 8)


DimPlot(tm_merged, pt.size = 0.5, reduction = 'umap', label = T, label.size = 2) + 
  scale_color_manual(values = cols[-c(6)])
ggsave('../tabulamuris_merged/umap/umap.seed210608156.1_2.celltype_refined.pdf', units = 'cm', width = 13.5, height = 8)
DimPlot(tm_merged, pt.size = 0.5, reduction = 'umap') + 
  scale_color_manual(values = cols[-c(6)]) + 
  theme_void() + 
  theme(legend.position = 'none')
ggsave('../tabulamuris_merged/umap/umap.seed210608156.1_2.celltype_refined.void.png', units = 'cm', width = 8, height = 8)




