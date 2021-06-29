library(Seurat)
library(ggplot2)
library(plyr)
library(dplyr)
library(harmony)

dir.create('tmp')
dir.create('stats')
dir.create('tsne')
dir.create('tsne/compare')
dir.create('umap')
dir.create('umap/compare')

counts <- readRDS('HCL_Fig1_counts.PB_BM.Rds')
label <- readRDS('HCL_Fig1_labels.PB_BM.Rds')
head(label)
identical(rownames(label), colnames(counts))

### The standard workflow ###
human_hcl <- CreateSeuratObject(counts = counts, project = "HCL")
human_hcl <- AddMetaData(object = human_hcl, metadata = label)
human_hcl@meta.data$orig.ident <- 'HCL'
head(human_hcl@meta.data); nrow(human_hcl@meta.data) # 21568 cells
summary(human_hcl@meta.data$batch)
summary(human_hcl@meta.data$donor)
remove(counts)

plt <- ggplot(human_hcl@meta.data, aes(batch, nFeature_RNA)) + geom_jitter(size = 0.25) + 
  theme_bw() + #ylim(0, 500) +
  theme(text = element_text(size = 7), axis.text = element_text(size = 7), panel.grid = element_blank(), axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5));plt
ggsave('stats/2-2.stats.human_hcl.nFeature_RNA.byLibrary.pdf', units = 'cm', width = 6, height = 6)
plt <- ggplot(human_hcl@meta.data, aes(batch, nCount_RNA)) + geom_jitter(size = 0.25) + 
  theme_bw() + #ylim(0, 800) +
  theme(text = element_text(size = 7), axis.text = element_text(size = 7), panel.grid = element_blank(), axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5));plt
ggsave('stats/2-2.stats.human_hcl.nCount_RNA.byLibrary.pdf', units = 'cm', width = 6, height = 6)

head(human_hcl@meta.data)
human_hcl <- NormalizeData(human_hcl, normalization.method = "LogNormalize", scale.factor = 10000)
human_hcl <- FindVariableFeatures(object = human_hcl, selection.method = "vst", nfeatures = 2000)
human_hcl <- ScaleData(object = human_hcl, vars.to.regress = c('nCount_RNA'))
human_hcl <- RunPCA(object = human_hcl, npcs = 60)
human_hcl <- JackStraw(object = human_hcl, num.replicate = 100, dims = 60)
human_hcl <- ScoreJackStraw(object = human_hcl, dims = 1:60)
JackStrawPlot(human_hcl, dims = 1:60) # 33 PCs 
ggsave('stats/2-4.pca.JackStrawPlot.mitoCut.pdf', units = 'cm', width = 25, height = 12)
saveRDS(human_hcl, 'tmp/human_hcl.Rds')

human_hcl <- RunHarmony(human_hcl, group.by.vars = "donor")
saveRDS(human_hcl, 'tmp/human_hcl.harmony.Rds')


### UMAP ###
for (seed in c(210320150:210320169)){
  human_hcl <- RunUMAP(human_hcl, dims=1:33, reduction = "harmony", reduction.key='UMAP', n.components=3, min.dist=0.3, seed.use = seed)
  DimPlot(human_hcl, group.by = "celltype", reduction = "umap", dims = c(1,2), label = T, label.size = 2)
  ggsave(paste0(c('umap/compare/umap.1_2.celltype.seed', seed, '.png'), collapse = ''), units = 'cm', width = 18, height = 10)
  DimPlot(human_hcl, group.by = "celltype", reduction = "umap", dims = c(1,3), label = T, label.size = 2)
  ggsave(paste0(c('umap/compare/umap.1_3.celltype.seed', seed, '.png'), collapse = ''), units = 'cm', width = 18, height = 10)
  DimPlot(human_hcl, group.by = "celltype", reduction = "umap", dims = c(2,3), label = T, label.size = 2)
  ggsave(paste0(c('umap/compare/umap.2_3.celltype.seed', seed, '.png'), collapse = ''), units = 'cm', width = 18, height = 10)
} # 210320163

### tSNE ###
for (seed in c(210320100:210320109)){
  human_hcl <- RunTSNE(human_hcl, dims=1:33, reduction = "harmony", reduction.key='tSNE', dim.embed=3, seed.use = seed)
  DimPlot(human_hcl, group.by = "celltype", reduction = "tsne", label = T, label.size = 2)
  ggsave(paste0(c('tsne/compare/tsne.1_2.celltype.seed', seed, '.png'), collapse = ''), units = 'cm', width = 18, height = 10)
} # 210320101

### UMAP fix ###
head(human_hcl@meta.data)
human_hcl <- RunUMAP(human_hcl, dims=1:33, reduction = "harmony", reduction.key='UMAP', n.components=3, min.dist=0.3, seed.use = 210320163)
DimPlot(object=human_hcl, group.by = "celltype", reduction = 'umap')
ggsave('umap/umap.combined.1_2.celltypes.seed210320163.pdf', units = 'cm', width = 18, height = 10)
DimPlot(object=human_hcl, group.by = "donor", reduction = 'umap')
ggsave('umap/umap.combined.1_2.libraries.seed210320163.pdf', units = 'cm', width = 13, height = 10)
FeaturePlot(human_hcl, features = 'nFeature_RNA', cols = c("grey","red"), reduction = "umap")
ggsave('umap/umap.combined.1_2.nFeature_RNA.seed210320163.pdf', units = 'cm', width = 12, height = 10)
FeaturePlot(human_hcl, features = 'nCount_RNA', cols = c("grey","red"), reduction = "umap")
ggsave('umap/umap.combined.1_2.nCount_RNA.seed210320163.pdf', units = 'cm', width = 12, height = 10)

DimPlot(object=human_hcl, group.by = "celltype", reduction = 'umap', dims = c(1,3))
ggsave('umap/umap.combined.1_3.celltypes.seed210320163.pdf', units = 'cm', width = 18, height = 10)
DimPlot(object=human_hcl, group.by = "donor", reduction = 'umap', dims = c(1,3))
ggsave('umap/umap.combined.1_3.libraries.seed210320163.pdf', units = 'cm', width = 13, height = 10)
FeaturePlot(human_hcl, features = 'nFeature_RNA', cols = c("grey","red"), reduction = "umap", dims = c(1,3))
ggsave('umap/umap.combined.1_3.nFeature_RNA.seed210320163.pdf', units = 'cm', width = 12, height = 10)
FeaturePlot(human_hcl, features = 'nCount_RNA', cols = c("grey","red"), reduction = "umap", dims = c(1,3))
ggsave('umap/umap.combined.1_3.nCount_RNA.seed210320163.pdf', units = 'cm', width = 12, height = 10)

DimPlot(object=human_hcl, group.by = "celltype", reduction = 'umap', dims = c(2,3))
ggsave('umap/umap.combined.2_3.celltypes.seed210320163.pdf', units = 'cm', width = 18, height = 10)
DimPlot(object=human_hcl, group.by = "donor", reduction = 'umap', dims = c(2,3))
ggsave('umap/umap.combined.2_3.libraries.seed210320163.pdf', units = 'cm', width = 13, height = 10)
FeaturePlot(human_hcl, features = 'nFeature_RNA', cols = c("grey","red"), reduction = "umap", dims = c(2,3))
ggsave('umap/umap.combined.2_3.nFeature_RNA.seed210320163.pdf', units = 'cm', width = 12, height = 10)
FeaturePlot(human_hcl, features = 'nCount_RNA', cols = c("grey","red"), reduction = "umap", dims = c(2,3))
ggsave('umap/umap.combined.2_3.nCount_RNA.seed210320163.pdf', units = 'cm', width = 12, height = 10)


### tSNE fix ###
human_hcl <- RunTSNE(human_hcl, dims=1:33, reduction = "harmony", reduction.key='tSNE', dim.embed=3, seed.use = 210320101)
DimPlot(object=human_hcl, group.by = "celltype", reduction = 'tsne')
ggsave('tsne/tsne.combined.1_2.celltypes.seed210320101.pdf', units = 'cm', width = 18, height = 10)
DimPlot(object=human_hcl, group.by = "donor", reduction = 'tsne')
ggsave('tsne/tsne.combined.1_2.libraries.seed210320101.pdf', units = 'cm', width = 13, height = 10)
FeaturePlot(human_hcl, features = 'nFeature_RNA', cols = c("grey","red"), reduction = "tsne")
ggsave('tsne/tsne.combined.1_2.nFeature_RNA.seed210320101.pdf', units = 'cm', width = 12, height = 10)
FeaturePlot(human_hcl, features = 'nCount_RNA', cols = c("grey","red"), reduction = "tsne")
ggsave('tsne/tsne.combined.1_2.nCount_RNA.seed210320101.pdf', units = 'cm', width = 12, height = 10)

Idents(human_hcl) <- 'celltype'
saveRDS(human_hcl, 'tmp/human_hcl.harmony.Rds')

dir.create('degs')
markers <- FindAllMarkers(object = human_hcl, min.pct = 0.25, logfc.threshold = 0.25, only.pos = T)
markers %>% group_by(cluster) %>% top_n(2, avg_logFC)
write.table(markers, 'degs/markers.human_hcl.celltype.txt', sep = '\t', quote = F, col.names = T, row.names = F) 
top10 <- markers %>% group_by(cluster) %>% top_n(10, avg_logFC)
DoHeatmap(object = human_hcl, features = top10$gene, angle = 90, size = 5, raster = T, draw.lines = F)
ggsave('degs/markers.human_hcl.celltype.pdf', units = 'cm', width = 40, height = 30)



###
library(RColorBrewer)
head(human_hcl@meta.data)
summary(human_hcl@meta.data$sample)
summary(human_hcl@meta.data$donor)
summary(human_hcl@meta.data$celltype)

my_cols <- colorRampPalette(brewer.pal(9, "Spectral"))
levels(human_hcl@meta.data$celltype)
human_hcl@meta.data$celltype <- factor(human_hcl@meta.data$celltype, levels = c("CB CD34+", "Monocyte", "Macrophage", "Dendritic cell", "B cell", "B cell (Plasmocyte)", 
                                                                                "T cell", "Erythroid progenitor cell (RP high)", "Erythroid cell", "Neutrophil (RPS high)", "Neutrophil"))

DimPlot(object=human_hcl, group.by = "celltype", reduction = 'umap', label = T, label.size = 2) + scale_color_manual(values = my_cols(11))
ggsave('umap/umap.combined.1_2.celltypes.seed210320163.pdf', units = 'cm', width = 18, height = 10)
DimPlot(object=human_hcl, group.by = "sample", reduction = 'umap') + scale_color_manual(values = c('red2', 'forestgreen'))
ggsave('umap/umap.combined.1_2.sample.seed210320163.pdf', units = 'cm', width = 15, height = 10)
DimPlot(object=human_hcl, group.by = "donor", reduction = 'umap') + scale_color_manual(values = colorRampPalette(brewer.pal(9, "Reds"))(6))
ggsave('umap/umap.combined.1_2.libraries.seed210320163.pdf', units = 'cm', width = 13, height = 10)


human_hcl@meta.data$celltype_broad <- mapvalues(human_hcl@meta.data$celltype, from = levels(human_hcl@meta.data$celltype),
                                                to = c("Progenitors", "Monocytes", "Macrophage", "DC", "B cells", "Plasma cells", 
                                                       "T cells", "Erythroids", "Erythroids", "Neutrophils", "Neutrophils"))
human_hcl@meta.data$celltype_broad <- factor(human_hcl@meta.data$celltype_broad, 
                                             levels = c("Progenitors", "Neutrophils", "Erythroids", "Monocytes", "Macrophage", "DC", "T cells", "B cells", "Plasma cells"))

DimPlot(object=human_hcl, group.by = "celltype_broad", reduction = 'umap', label = T, label.size = 2) + scale_color_manual(values = my_cols(9))
ggsave('umap/umap.combined.1_2.celltype_broad.seed210320163.pdf', units = 'cm', width = 13.5, height = 10)
Idents(human_hcl) <- 'celltype_broad'
#saveRDS(human_hcl, 'tmp/human_hcl.harmony.Rds')

markers <- FindAllMarkers(object = human_hcl, min.pct = 0.25, logfc.threshold = 0.25, only.pos = T)
markers %>% group_by(cluster) %>% top_n(2, avg_logFC)
write.table(markers, 'degs/markers.human_hcl.celltype_broad.txt', sep = '\t', quote = F, col.names = T, row.names = F) 
top10 <- markers %>% group_by(cluster) %>% top_n(10, avg_logFC)
DoHeatmap(object = human_hcl, features = top10$gene, angle = 90, size = 5, raster = T, draw.lines = F)
ggsave('degs/markers.human_hcl.celltype_broad.pdf', units = 'cm', width = 40, height = 30)



### Clustering - refined ###
dir.create('res')
dir.create('res/compare')
human_hcl <- FindNeighbors(object = human_hcl, dims = 1:33)
for (res in c(1:12)){
  res <- res/10
  human_hcl <- FindClusters(object = human_hcl, resolution = res)
  
  DimPlot(human_hcl, label = T, pt.size = 0.5, reduction = 'umap')
  ggsave(paste(c('res/compare/umap.seed210320163.1_2.res_',res,'.pdf'), collapse = ''), units = 'cm', width = 12.5, height = 10)
  DimPlot(human_hcl, label = T, pt.size = 0.5, reduction = 'tsne')
  ggsave(paste(c('res/compare/tsne.seed210320101.1_2.res_',res,'.pdf'), collapse = ''), units = 'cm', width = 12.5, height = 10)
}

human_hcl <- readRDS('tmp/human_hcl.harmony.Rds')
human_hcl <- FindNeighbors(object = human_hcl, dims = 1:33)

human_hcl <- FindClusters(object = human_hcl, resolution = 0.3)
DimPlot(human_hcl, label = T, pt.size = 0.5, reduction = 'umap')
ggsave('res/umap.seed210320163.1_2.res_0.3.pdf', units = 'cm', width = 11.5, height = 10)
DimPlot(human_hcl, label = T, dims = c(1,3), pt.size = 0.5, reduction = 'umap')
ggsave('res/umap.seed210320163.1_3.res_0.3.pdf', units = 'cm', width = 11.5, height = 10)
DimPlot(human_hcl, label = T, dims = c(2,3), pt.size = 0.5, reduction = 'umap')
ggsave('res/umap.seed210320163.2_3.res_0.3.pdf', units = 'cm', width = 11.5, height = 10)

DimPlot(human_hcl, label = T, pt.size = 0.5, reduction = 'tsne')
ggsave('res/tsne.seed210320101.1_2.res_0.3.pdf', units = 'cm', width = 11.5, height = 10)
DimPlot(human_hcl, label = T, dims = c(1,3), pt.size = 0.5, reduction = 'tsne')
ggsave('res/tsne.seed210320101.1_3.res_0.3.pdf', units = 'cm', width = 11.5, height = 10)
DimPlot(human_hcl, label = T, dims = c(2,3), pt.size = 0.5, reduction = 'tsne')
ggsave('res/tsne.seed210320101.2_3.res_0.3.pdf', units = 'cm', width = 11.5, height = 10)
saveRDS(human_hcl, 'tmp/human_hcl.harmony.Rds')

markers <- FindAllMarkers(object = human_hcl, min.pct = 0.25, logfc.threshold = 0.25, only.pos = T)
markers %>% group_by(cluster) %>% top_n(2, avg_logFC)
write.table(markers, 'degs/markers.human_hcl.res_0.3.txt', sep = '\t', quote = F, col.names = T, row.names = F) 
top10 <- markers %>% group_by(cluster) %>% top_n(10, avg_logFC)
DoHeatmap(object = human_hcl, features = top10$gene, angle = 90, size = 5, raster = T, draw.lines = F)
ggsave('degs/markers.human_hcl.res_0.3.pdf', units = 'cm', width = 40, height = 30)


### 
head(human_hcl@meta.data)
DimPlot(human_hcl) + facet_wrap(~human_hcl@meta.data$RNA_snn_res.0.3)
human_hcl@meta.data$seurat_clusters <- NULL
levels(human_hcl@meta.data$RNA_snn_res.0.3)

VlnPlot(human_hcl, features = c('CD3D', 'CD8A', 'CD4', 'NKG7', 'KLRD1'), pt.size = 0)
VlnPlot(human_hcl, features = c('FCGR3B', 'ALPL', 'CXCR1', 'CXCR2', 'ADGRG3', 'CMTM2', 'PROK2', 'MME', 'MMP25'), pt.size = 0)
VlnPlot(human_hcl, features = c('CD34', 'PRSS57', 'MPO', 'AZU1', 'ELANE', 'PRTN3'), pt.size = 0)
VlnPlot(human_hcl, features = c('S100A9', 'S100A8', 'S100A12', 'VCAN', 'FCN1', 'S100A6'), pt.size = 0)
#  0  1  2  3  4
#  5  6  7  8  9
# 10 11 12 13

celltypes <- c('Monocytes', 'Erythroids', 'T cells', 'Granulocyte progenitor', 'NK cells',
               'Neutrophils', 'Neutrophils', 'Macrophage', 'DC', 'Plasma cell',
               'T cells', 'Monocytes', 'B cells', 'pDC')


human_hcl@meta.data$celltype_refined <- mapvalues(human_hcl@meta.data$RNA_snn_res.0.3, from = c(0:13), to = celltypes)
human_hcl@meta.data$celltype_refined <- factor(human_hcl@meta.data$celltype_refined, 
                                               levels = c("Granulocyte progenitor", "Neutrophils", "Erythroids", "Monocytes", "Macrophage", "DC", "pDC", 
                                                          "T cells", "NK cells", "B cells", "Plasma cell"))
levels(human_hcl@meta.data$celltype_refined)
Idents(human_hcl) <- 'celltype_refined'
saveRDS(human_hcl, 'tmp/human_hcl.harmony.Rds')


DimPlot(human_hcl, label = T, pt.size = 0.5, reduction = 'umap')
ggsave('res/umap.seed210320163.1_2.celltype_refined.pdf', units = 'cm', width = 15, height = 10)
DimPlot(human_hcl, label = T, dims = c(1,3), pt.size = 0.5, reduction = 'umap')
ggsave('res/umap.seed210320163.1_3.celltype_refined.pdf', units = 'cm', width = 15, height = 10)
DimPlot(human_hcl, label = T, dims = c(2,3), pt.size = 0.5, reduction = 'umap')
ggsave('res/umap.seed210320163.2_3.celltype_refined.pdf', units = 'cm', width = 15, height = 10)

DimPlot(human_hcl, label = T, pt.size = 0.5, reduction = 'tsne')
ggsave('res/tsne.seed210320101.1_2.celltype_refined.pdf', units = 'cm', width = 15, height = 10)
DimPlot(human_hcl, label = T, dims = c(1,3), pt.size = 0.5, reduction = 'tsne')
ggsave('res/tsne.seed210320101.1_3.celltype_refined.pdf', units = 'cm', width = 15, height = 10)
DimPlot(human_hcl, label = T, dims = c(2,3), pt.size = 0.5, reduction = 'tsne')
ggsave('res/tsne.seed210320101.2_3.celltype_refined.pdf', units = 'cm', width = 15, height = 10)


markers <- FindAllMarkers(object = human_hcl, min.pct = 0.25, logfc.threshold = 0.25, only.pos = T)
markers %>% group_by(cluster) %>% top_n(2, avg_logFC)
write.table(markers, 'degs/markers.human_hcl.celltype_refined.txt', sep = '\t', quote = F, col.names = T, row.names = F) 
top10 <- markers %>% group_by(cluster) %>% top_n(10, avg_logFC)
DoHeatmap(object = human_hcl, features = top10$gene, angle = 90, size = 5, raster = T, draw.lines = F)
ggsave('degs/markers.human_hcl.celltype_refined.pdf', units = 'cm', width = 40, height = 30)




library(Seurat)
library(plyr)
human_hcl <- readRDS('tmp/human_hcl.harmony.Rds')
head(human_hcl@meta.data); nrow(human_hcl@meta.data)

### Pseudo-cell (10 cells) transformation ###
pseudocell_total <- data.frame(matrix(ncol = 0, nrow = length(rownames(human_hcl))))
rownames(pseudocell_total) <- rownames(human_hcl)

for (ct in levels(human_hcl@meta.data$celltype_refined)) {
  tmpObj <- subset(human_hcl, celltype_refined == ct)
  tmpObj_exprs <- as.matrix(GetAssayData(tmpObj, slot = 'data', assay = 'RNA'))
  
  bcs <- rownames(tmpObj@meta.data)
  rounds <- trunc(length(bcs)/10)
  
  ct_pseudocells <- data.frame(matrix(ncol = rounds, nrow = length(rownames(human_hcl))))
  colnames(ct_pseudocells) <- unlist(lapply(c(1:rounds), function (x) paste0(c(as.character(ct), x), collapse = ' pc') ))
  rownames(ct_pseudocells) <- rownames(human_hcl)
  
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


###
library(Seurat)
library(plyr)
library(ggplot2)
library(RColorBrewer)

human_hcl <- readRDS('tmp/human_hcl.harmony.Rds')
human_hcl@meta.data$celltype_refined <- factor(human_hcl@meta.data$celltype_refined, 
                                               levels = c('Granulocyte progenitor', 'Erythroids', 'Neutrophils', 'Monocytes', 
                                                          'Macrophage', 'DC', 'pDC', 'T cells', 'NK cells', 'B cells', 'Plasma cell'))
Idents(human_hcl) <- 'celltype_refined'
head(human_hcl@meta.data); nrow(human_hcl@meta.data)

mycol <- colorRampPalette(brewer.pal(9, "Spectral"))(17)
mycol
# hcl mycol[-c(1, 10, 12, 13, 14, 17)]
#   1       2           3     4       5       6     7     8       9         10    11        12        13        14        15        16    17
#         gran.prog    Ery    Neu    Mono    Mac    DC    pDC    T cell           NK cell                                B cell    PC

# hca mycol[-c(4, 6)]
#   1       2           3     4       5       6     7     8       9         10    11        12        13        14        15        16    17
# Prog    gran.prog    Ery           Mono           DC    pDC    T cell    CD8    NK cell    pre-PC    pro-B    pre-B    B cell    PC    platelet


DimPlot(human_hcl, pt.size = 0.5, reduction = 'umap') + 
  scale_color_manual(values =  mycol[-c(1, 10, 12, 13, 14, 17)] ) + 
  theme_void() + 
  theme(legend.position = 'none')
ggsave('res/umap.seed210320163.1_2.celltype_refined.void.png', units = 'cm', width = 8, height = 8)

DimPlot(human_hcl, pt.size = 0.5, reduction = 'umap', label = T, label.size = 2) + 
  scale_color_manual(values = mycol[-c(1, 10, 12, 13, 14, 17)] )
ggsave('res/umap.seed210320163.1_2.celltype_refined.pdf', units = 'cm', width = 13.5, height = 8)


