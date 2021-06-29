library(Seurat)
library(ggplot2)
library(plyr)
library(dplyr)
library(harmony)

label <- readRDS('../merged/label.Rds')
head(label)

### tabula muris 10X ###
tm_tenx <- readRDS('/home/sangho/2020_newbuild/projects/drosophila/projectX_review/Cross-species/mouse-TabulaMuris/droplet_Marrow_seurat_tiss.counts.Rds')
tm_tenx_bcs <- intersect(rownames(subset(label, dataset == '10X')), colnames(tm_tenx))

tm_tenx <- tm_tenx[, tm_tenx_bcs]
tm_tenx$sums <- rowSums(tm_tenx)
tm_tenx <- subset(tm_tenx, sums != 0)
tm_tenx$sums <- NULL
mouse <- CreateSeuratObject(counts = tm_tenx, project = "Tabula Muris"); mouse@meta.data$orig.ident <- 'Tabula Muris'
mouse <- AddMetaData(object = mouse, metadata = label[tm_tenx_bcs, ])
mouse@meta.data <- droplevels(mouse@meta.data)
head(mouse@meta.data); summary(mouse@meta.data)

dir.create('tmp')
saveRDS(mouse, 'tmp/tm_10x.Rds')

### Standard workflow ###
dir.create('stats')
dir.create('tsne')
dir.create('tsne/compare')
dir.create('umap')
dir.create('umap/compare')
dir.create('degs')

#mouse <- readRDS('tmp/tm_10x.Rds')
head(mouse@meta.data); nrow(mouse@meta.data) # 3652 cells

mouse <- NormalizeData(mouse, normalization.method = "LogNormalize", scale.factor = 10000)
mouse <- FindVariableFeatures(object = mouse, selection.method = "vst", nfeatures = 2000)
mouse <- ScaleData(object = mouse, vars.to.regress = c('nCount_RNA'))
mouse <- RunPCA(object = mouse, npcs = 60)
mouse <- JackStraw(object = mouse, num.replicate = 100, dims = 60)
mouse <- ScoreJackStraw(object = mouse, dims = 1:60)
JackStrawPlot(mouse, dims = 1:60) # 38 PCs 
ggsave('stats/2-4.pca.JackStrawPlot.mitoCut.pdf', units = 'cm', width = 25, height = 12)
saveRDS(mouse, 'tmp/tm_10x_2.Rds')

mouse <- RunHarmony(mouse, group.by.vars = "mouse.id")
saveRDS(mouse, 'tmp/tm_10x.harmony.Rds')

### UMAP ###
for (seed in c(210517250:210517269)){
  mouse <- RunUMAP(mouse, dims=1:38, reduction = "harmony", reduction.key='UMAP', n.components=3, min.dist=0.3, seed.use = seed)
  DimPlot(mouse, group.by = "celltype_merged", reduction = "umap", label = T, label.size = 2)
  ggsave(paste0(c('umap/compare/umap.1_2.markers.seed', seed, '.png'), collapse = ''), units = 'cm', width = 14, height = 10)
} # 210517250

### tSNE ###
for (seed in c(210517200:210517209)){
  mouse <- RunTSNE(mouse, dims=1:38, reduction = "harmony", reduction.key='tSNE', dim.embed=3, seed.use = seed)
  DimPlot(mouse, group.by = "celltype_merged", reduction = "tsne", label = T, label.size = 2)
  ggsave(paste0(c('tsne/compare/tsne.1_2.markers.seed', seed, '.png'), collapse = ''), units = 'cm', width = 14, height = 10)
} # 210517201


mouse <- RunUMAP(mouse, dims=1:38, reduction = "harmony", reduction.key='UMAP', n.components=3, min.dist=0.3, seed.use = 210517250)
DimPlot(object=mouse, group.by = "dataset", reduction = 'umap')
ggsave('umap/umap.1_2.dataset.seed210517250.pdf', units = 'cm', width = 13.5, height = 10)
DimPlot(object=mouse, group.by = "tissue", reduction = 'umap')
ggsave('umap/umap.1_2.tissue.seed210517250.pdf', units = 'cm', width = 14.5, height = 10)
DimPlot(object=mouse, group.by = "celltype_merged", reduction = 'umap', label = T, label.size = 2)
ggsave('umap/umap.1_2.celltype_merged.seed210517250.pdf', units = 'cm', width = 14.5, height = 10)

mouse <- RunTSNE(mouse, dims=1:38, reduction = "harmony", reduction.key='tSNE', dim.embed=3, seed.use = 210517201)
DimPlot(object=mouse, group.by = "dataset", reduction = 'tsne')
ggsave('tsne/tsne.1_2.dataset.seed210517201.pdf', units = 'cm', width = 13.5, height = 10)
DimPlot(object=mouse, group.by = "tissue", reduction = 'tsne')
ggsave('tsne/tsne.1_2.tissue.seed210517201.pdf', units = 'cm', width = 14.5, height = 10)
DimPlot(object=mouse, group.by = "celltype_merged", reduction = 'tsne', label = T, label.size = 2)
ggsave('tsne/tsne.1_2.celltype_merged.seed210517201.pdf', units = 'cm', width = 14.5, height = 10)
saveRDS(mouse, 'tmp/tm_10x.harmony.Rds')


### clustering ###
dir.create('res')
dir.create('res/compare')
mouse <- FindNeighbors(object = mouse, dims = 1:38)
for (res in c(1:12)){
  res <- res/10
  mouse <- FindClusters(object = mouse, resolution = res)
  DimPlot(mouse, label = T, pt.size = 0.5, reduction = 'umap')
  ggsave(paste(c('res/compare/umap.seed210517250.1_2.res_',res,'.pdf'), collapse = ''), units = 'cm', width = 13.5, height = 10)
  DimPlot(mouse, label = T, dims = c(1,3), pt.size = 0.5, reduction = 'umap')
  ggsave(paste(c('res/compare/umap.seed210517250.1_3.res_',res,'.pdf'), collapse = ''), units = 'cm', width = 13.5, height = 10)
  DimPlot(mouse, label = T, dims = c(2,3), pt.size = 0.5, reduction = 'umap')
  ggsave(paste(c('res/compare/umap.seed210517250.2_3.res_',res,'.pdf'), collapse = ''), units = 'cm', width = 13.5, height = 10)
  
  DimPlot(mouse, label = T, pt.size = 0.5, reduction = 'tsne')
  ggsave(paste(c('res/compare/tsne.seed210517201.1_2.res_',res,'.pdf'), collapse = ''), units = 'cm', width = 13.5, height = 10)
  DimPlot(mouse, label = T, dims = c(1,3), pt.size = 0.5, reduction = 'tsne')
  ggsave(paste(c('res/compare/tsne.seed210517201.1_3.res_',res,'.pdf'), collapse = ''), units = 'cm', width = 13.5, height = 10)
  DimPlot(mouse, label = T, dims = c(2,3), pt.size = 0.5, reduction = 'tsne')
  ggsave(paste(c('res/compare/tsne.seed210517201.2_3.res_',res,'.pdf'), collapse = ''), units = 'cm', width = 13.5, height = 10)
}

mouse <- readRDS('tmp/tm_10x.harmony.Rds')
mouse <- FindNeighbors(object = mouse, dims = 1:38)

mouse <- FindClusters(object = mouse, resolution = 0.5)
DimPlot(mouse, label = T, pt.size = 0.5, reduction = 'umap')
ggsave('res/umap.seed210517250.1_2.res_0.5.pdf', units = 'cm', width = 12.5, height = 10)
DimPlot(mouse, label = T, dims = c(1,3), pt.size = 0.5, reduction = 'umap')
ggsave('res/umap.seed210517250.1_3.res_0.5.pdf', units = 'cm', width = 12.5, height = 10)
DimPlot(mouse, label = T, dims = c(2,3), pt.size = 0.5, reduction = 'umap')
ggsave('res/umap.seed210517250.2_3.res_0.5.pdf', units = 'cm', width = 12.5, height = 10)

DimPlot(mouse, label = T, pt.size = 0.5, reduction = 'tsne')
ggsave('res/tsne.seed210517201.1_2.res_0.5.pdf', units = 'cm', width = 12.5, height = 10)
DimPlot(mouse, label = T, dims = c(1,3), pt.size = 0.5, reduction = 'tsne')
ggsave('res/tsne.seed210517201.1_3.res_0.5.pdf', units = 'cm', width = 12.5, height = 10)
DimPlot(mouse, label = T, dims = c(2,3), pt.size = 0.5, reduction = 'tsne')
ggsave('res/tsne.seed210517201.2_3.res_0.5.pdf', units = 'cm', width = 12.5, height = 10)
saveRDS(mouse, 'tmp/tm_10x.harmony.Rds')


### markers ###
markers <- FindAllMarkers(object = mouse, min.pct = 0.25, logfc.threshold = 0.25, only.pos = T)
markers %>% group_by(cluster) %>% top_n(2, avg_logFC)
write.table(markers, 'degs/markers.mouse.res_0.5.txt', sep = '\t', quote = F, col.names = T, row.names = F) 
top10 <- markers %>% group_by(cluster) %>% top_n(10, avg_logFC)
DoHeatmap(object = mouse, features = top10$gene, angle = 90, size = 5, raster = T, draw.lines = F)
ggsave('degs/markers.mouse.res_0.5.pdf', units = 'cm', width = 30, height = 20)




library(Seurat)
library(ggplot2)
library(plyr)
library(dplyr)

mouse <- readRDS('tmp/tm_10x.harmony.Rds')
mouse <- readRDS('tmp/tm_10x.harmony.clustered.Rds')

DimPlot(mouse, label = T)

VlnPlot(mouse, pt.size = 0, features = c('S100a9', 'S100a8', 'Mmp9', 'Mmp8')) # Neutrophils
VlnPlot(mouse, pt.size = 0, features = c('Gzma', 'Klrb1c', 'Ncr1')) # NK cells
VlnPlot(mouse, pt.size = 0, features = c('Cd8a', 'Cd3e')) # T cells
VlnPlot(mouse, pt.size = 0, features = c('Cd79a', 'Cd79b', 'Mzb1', 'Jchain', 'Iglv1', 'Iglc1', 'Iglc2', 'Derl3', 'Tnfrsf17')) # B PC
VlnPlot(mouse, pt.size = 0, features = c('Prss34', 'Mcpt8', 'Ccl3', 'Cpa3', 'Ccl4', 'Fcer1a')) # Mast cells
VlnPlot(mouse, pt.size = 0, features = c('Siglech', 'Ccr9', 'Bst2', 'Pacsin1', 'Tcf4')) # pDC
VlnPlot(mouse, pt.size = 0, features = c('Csf1r', 'C1qa', 'Apoe', 'Fn1')) # Mono macrophage
VlnPlot(mouse, pt.size = 0, features = c('Alas1', 'Anxa3', 'Elane', 'Gatm', 'Mpo', 'Hk3', 'Hp', 'Igsf6')) # GMP
VlnPlot(mouse, pt.size = 0, features = c('Kit', 'Cd34')) # 


levels(mouse@meta.data$RNA_snn_res.0.5)
### rename ### 
# 0 1 2 3 4 5
# 6 7 8 9 10 11
# 12 13 14 15

clusters <- c('Neutrophils', 'Monocytes', 'Granulocytopoietic cells', 'Progenitors', 'B cells', 'Erythroids',
              'Monocytes', 'pDC', 'T cells', 'Erythroids', 'Neutrophils', 'B cells',
              'DC', 'B cells', 'Basophil', 'Plasma cells')
mouse@meta.data$celltype_refined <- mapvalues(mouse@meta.data$RNA_snn_res.0.5, from = levels(mouse@meta.data$RNA_snn_res.0.5), to = clusters)
mouse@meta.data$celltype_refined <- factor(mouse@meta.data$celltype_refined, 
                                           levels = c("Progenitors", "Granulocytopoietic cells", "Neutrophils", "Erythroids", "Monocytes", "DC", "pDC", "Basophil", "T cells", "B cells", "Plasma cells"))
Idents(mouse) <- 'celltype_refined'
saveRDS(mouse, 'tmp/tm_10x.harmony.clustered.Rds')


DimPlot(mouse, label = T, pt.size = 0.5, reduction = 'umap', label.size = 2)
ggsave('res/umap.seed210517250.1_2.celltype_refined.pdf', units = 'cm', width = 16, height = 10)
DimPlot(mouse, label = T, dims = c(1,3), pt.size = 0.5, reduction = 'umap', label.size = 2)
ggsave('res/umap.seed210517250.1_3.celltype_refined.pdf', units = 'cm', width = 16, height = 10)
DimPlot(mouse, label = T, dims = c(2,3), pt.size = 0.5, reduction = 'umap', label.size = 2)
ggsave('res/umap.seed210517250.2_3.celltype_refined.pdf', units = 'cm', width = 16, height = 10)

DimPlot(mouse, label = T, pt.size = 0.5, reduction = 'tsne', label.size = 2)
ggsave('res/tsne.seed210517201.1_2.celltype_refined.pdf', units = 'cm', width = 16, height = 10)
DimPlot(mouse, label = T, dims = c(1,3), pt.size = 0.5, reduction = 'tsne', label.size = 2)
ggsave('res/tsne.seed210517201.1_3.celltype_refined.pdf', units = 'cm', width = 16, height = 10)
DimPlot(mouse, label = T, dims = c(2,3), pt.size = 0.5, reduction = 'tsne', label.size = 2)
ggsave('res/tsne.seed210517201.2_3.celltype_refined.pdf', units = 'cm', width = 16, height = 10)


library(Seurat)
library(plyr)
library(dplyr)
library(ggplot2)

mouse <- readRDS('tmp/tm_10x.harmony.clustered.Rds')
levels(Idents(mouse))

markers <- FindAllMarkers(object = mouse, min.pct = 0.25, logfc.threshold = 0.25, only.pos = T, test.use = 'MAST')
markers %>% group_by(cluster) %>% top_n(2, avg_logFC)
write.table(markers, 'degs/markers.mouse.10x.celltype_refined.MAST.txt', sep = '\t', quote = F, col.names = T, row.names = F) 
top10 <- markers %>% group_by(cluster) %>% top_n(10, avg_logFC)
DoHeatmap(object = mouse, features = top10$gene, angle = 90, size = 3, raster = T, draw.lines = F)
ggsave('degs/markers.mouse.10x.celltype_refined.MAST.pdf', units = 'cm', width = 40, height = 30)



library(Seurat)
library(plyr)
library(dplyr)
library(ggplot2)
mouse <- readRDS('tmp/tm_10x.harmony.clustered.Rds')
head(mouse@meta.data)

### Pseudo-cell (10 cells) transformation ###
pseudocell_total <- data.frame(matrix(ncol = 0, nrow = length(rownames(mouse))))
rownames(pseudocell_total) <- rownames(mouse)

for (ct in levels(mouse@meta.data$celltype_refined)) {
  tmpObj <- subset(mouse, celltype_refined == ct)
  tmpObj_exprs <- as.matrix(GetAssayData(tmpObj, slot = 'data', assay = 'RNA'))
  
  bcs <- rownames(tmpObj@meta.data)
  rounds <- trunc(length(bcs)/10)
  
  ct_pseudocells <- data.frame(matrix(ncol = rounds, nrow = length(rownames(mouse))))
  colnames(ct_pseudocells) <- unlist(lapply(c(1:rounds), function (x) paste0(c(as.character(ct), x), collapse = ' pc') ))
  rownames(ct_pseudocells) <- rownames(mouse)
  
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
saveRDS(pseudocell_total, 'tmp/celltype_refined.pseudocell10.Rds')


###
library(Seurat)
library(plyr)
library(dplyr)
library(ggplot2)
library(RColorBrewer)

mouse <- readRDS('tmp/tm_10x.harmony.clustered.Rds')
head(mouse@meta.data)

DimPlot(mouse, pt.size = 0.5, reduction = 'umap') + 
  scale_color_manual(values = colorRampPalette(brewer.pal(9, "Spectral"))(length(levels(mouse@meta.data$celltype_refined))) ) + 
  theme_void() + 
  theme(legend.position = 'none')
ggsave('res/umap.seed210517250.1_2.celltype_refined.void.png', units = 'cm', width = 8, height = 8)
DimPlot(mouse, pt.size = 0.5, reduction = 'umap', label = T, label.size = 2) + 
  scale_color_manual(values = colorRampPalette(brewer.pal(9, "Spectral"))(length(levels(mouse@meta.data$celltype_refined))) )
ggsave('res/umap.seed210517250.1_2.celltype_refined.pdf', units = 'cm', width = 13.5, height = 8)


