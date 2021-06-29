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

exprs <- 'GSE100910_indrop_counts.txt'

data <- read.delim(exprs, check.names = F, row.names = 1)
data[1:4, 1:4]
rownames(data) <- unlist(lapply(rownames(data), function (x) unlist(strsplit(as.character(x), split = '_1_of_many'))[1] ))

labelData <- data.frame(row.names = colnames(data), Library = unlist(lapply(colnames(data), function (x) unlist(strsplit(as.character(x), split = ':'))[2] )))
levels(labelData$Library); summary(labelData$Library)
labelData <- subset(labelData, Library %in% c('WT_1', 'WT_2', 'WT_3'))
labelData <- droplevels(labelData)
levels(labelData$Library)
head(labelData)

data <- data[, rownames(labelData)]; dim(data)

dropseq <- CreateSeuratObject(counts = data, project = "zebrafish_InDrop")
dropseq <- AddMetaData(object = dropseq, metadata = labelData)
head(dropseq@meta.data); nrow(dropseq@meta.data) # 3782 cells

plt <- ggplot(dropseq@meta.data, aes(Library, nFeature_RNA)) + geom_jitter(size = 0.25) + 
  theme_bw() + 
  theme(text = element_text(size = 7), axis.text = element_text(size = 7), panel.grid = element_blank(), axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5));plt
ggsave('stats/2-2.stats.nFeature_RNA.byLibrary.pdf', units = 'cm', width = 6, height = 8)
plt <- ggplot(dropseq@meta.data, aes(Library, nCount_RNA)) + geom_jitter(size = 0.25) + 
  theme_bw() + 
  theme(text = element_text(size = 7), axis.text = element_text(size = 7), panel.grid = element_blank(), axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5));plt
ggsave('stats/2-2.stats.nCount_RNA.byLibrary.pdf', units = 'cm', width = 6, height = 8)

### The standard workflow ###
dropseq.combined <- NormalizeData(dropseq, normalization.method = "LogNormalize", scale.factor = 10000)
dropseq.combined <- FindVariableFeatures(object = dropseq.combined, selection.method = "vst", nfeatures = 2000)

dropseq.combined <- ScaleData(object = dropseq.combined, vars.to.regress = c('Library', 'nCount_RNA'))
dropseq.combined <- RunPCA(object = dropseq.combined, npcs = 60)
dropseq.combined <- JackStraw(object = dropseq.combined, num.replicate = 100, dims = 60)
dropseq.combined <- ScoreJackStraw(object = dropseq.combined, dims = 1:60)
JackStrawPlot(dropseq.combined, dims = 1:60) #  PCs 
ggsave('stats/2-4.pca.JackStrawPlot.mitoCut.pdf', units = 'cm', width = 25, height = 12)
#saveRDS(dropseq.combined, 'tmp/zebrafish.Rds')

dropseq.combined <- RunHarmony(dropseq.combined, group.by.vars = "Library")
#saveRDS(dropseq.combined, 'tmp/zebrafish.harmony.Rds')


genes <- c('itga2b', 'rag1', 'epor', 'lyz', 'marco', 'grna', 'cd79a', 'cd8a', 'zap70'); genes %in% rownames(dropseq.combined)
### UMAP ###
for (seed in c(210218150:210218169)){
  dropseq.combined <- RunUMAP(dropseq.combined, dims=1:28, reduction = "harmony", reduction.key='UMAP', n.components=3, min.dist=0.3, seed.use = seed)
  FeaturePlot(dropseq.combined, features = genes, cols = c("grey","red"), reduction = "umap")
  ggsave(paste0(c('umap/compare/umap.1_2.markers.seed', seed, '.png'), collapse = ''), units = 'cm', width = 35, height = 30)
} # 210218150

### tSNE ###
for (seed in c(210218100:210218119)){
  dropseq.combined <- RunTSNE(dropseq.combined, dims=1:28, reduction = "harmony", reduction.key='tSNE', dim.embed=3, seed.use = seed)
  FeaturePlot(dropseq.combined, features = genes, cols = c("grey","red"), reduction = "tsne")
  ggsave(paste0(c('tsne/compare/tsne.1_2.markers.seed', seed, '.png'), collapse = ''), units = 'cm', width = 35, height = 30)
} # 210218101, 210218102

### UMAP fix ###
head(dropseq.combined@meta.data)
dropseq.combined <- RunUMAP(dropseq.combined, dims=1:28, reduction = "harmony", reduction.key='UMAP', n.components=3, min.dist=0.3, seed.use = 210218150)
FeaturePlot(dropseq.combined, features = genes, cols = c("grey","red"), reduction = "umap")
ggsave('umap/umap.combined.1_2.markers.seed210218150.png', units = 'cm', width = 35, height = 30)
DimPlot(object=dropseq.combined, group.by = "Library", reduction = 'umap')
ggsave('umap/umap.combined.1_2.libraries.seed210218150.pdf', units = 'cm', width = 12, height = 10)
FeaturePlot(dropseq.combined, features = 'nFeature_RNA', cols = c("grey","red"), reduction = "umap")
ggsave('umap/umap.combined.1_2.nFeature_RNA.seed210218150.pdf', units = 'cm', width = 12, height = 10)
FeaturePlot(dropseq.combined, features = 'nCount_RNA', cols = c("grey","red"), reduction = "umap")
ggsave('umap/umap.combined.1_2.nCount_RNA.seed210218150.pdf', units = 'cm', width = 12, height = 10)

### tSNE fix ###
dropseq.combined <- RunTSNE(dropseq.combined, dims=1:28, reduction = "harmony", reduction.key='tSNE', dim.embed=3, seed.use = 210218102)
FeaturePlot(dropseq.combined, features = genes, cols = c("grey","red"), reduction = "tsne")
ggsave('tsne/tsne.combined.1_2.markers.seed210218102.png', units = 'cm', width = 35, height = 30)
DimPlot(object=dropseq.combined, group.by = "Library", reduction = 'tsne')
ggsave('tsne/tsne.combined.1_2.libraries.seed210218102.pdf', units = 'cm', width = 12, height = 10)
FeaturePlot(dropseq.combined, features = 'nFeature_RNA', cols = c("grey","red"), reduction = "tsne")
ggsave('tsne/tsne.combined.1_2.nFeature_RNA.seed210218102.pdf', units = 'cm', width = 12, height = 10)
FeaturePlot(dropseq.combined, features = 'nCount_RNA', cols = c("grey","red"), reduction = "tsne")
ggsave('tsne/tsne.combined.1_2.nCount_RNA.seed210218102.pdf', units = 'cm', width = 12, height = 10)

#saveRDS(dropseq.combined, 'tmp/zebrafish.harmony.Rds')

### Clustering ###
#dir.create('res/')
#dir.create('res/compare/')

#dropseq.combined <- readRDS('tmp/zebrafish.harmony.Rds')
dropseq.combined <- FindNeighbors(object = dropseq.combined, dims = 1:28)
for (res in c(1:15)){
  res <- res/10
  dropseq.combined <- FindClusters(object = dropseq.combined, resolution = res)
  
  DimPlot(dropseq.combined, label = T, dims = c(1,2), reduction = 'tsne')
  ggsave(paste(c('res/compare/tsne.1_2.res_',res,'.png'), collapse = ''), units = 'cm', width = 12, height = 10)
  DimPlot(dropseq.combined, label = T, dims = c(1,2), reduction = 'umap')
  ggsave(paste(c('res/compare/umap.1_2.res_',res,'.png'), collapse = ''), units = 'cm', width = 12, height = 10)
}

dropseq.combined <- readRDS('tmp/zebrafish.harmony.Rds')
dropseq.combined <- FindNeighbors(object = dropseq.combined, dims = 1:28)
dropseq.combined <- FindClusters(object = dropseq.combined, resolution = 0.7)

DimPlot(dropseq.combined, label = T, pt.size = 0.5, reduction = 'tsne')
ggsave('res/tsne.1_2.res_0.7.pdf', units = 'cm', width = 12, height = 8.5)
DimPlot(dropseq.combined, label = T, dims = c(1,3), pt.size = 0.5, reduction = 'tsne')
ggsave('res/tsne.1_3.res_0.7.pdf', units = 'cm', width = 12, height = 8.5)
DimPlot(dropseq.combined, label = T, dims = c(2,3), pt.size = 0.5, reduction = 'tsne')
ggsave('res/tsne.2_3.res_0.7.pdf', units = 'cm', width = 12, height = 8.5)

DimPlot(dropseq.combined, reduction = "umap", pt.size = .5, label = T)
ggsave('res/umap.1_2.min_dist_0.3.res_0.7.pdf', units = 'cm', width = 12, height = 8.5)
DimPlot(dropseq.combined, reduction = "umap", pt.size = .5, dims = c(1, 3), label = T)
ggsave('res/umap.1_3.min_dist_0.3.res_0.7.pdf', units = 'cm', width = 12, height = 8.5)
DimPlot(dropseq.combined, reduction = "umap", pt.size = .5, dims = c(2, 3), label = T)
ggsave('res/umap.2_3.min_dist_0.3.res_0.7.pdf', units = 'cm', width = 12, height = 8.5)
#saveRDS(dropseq.combined, 'tmp/zebrafish.harmony.Rds')


### DEG analysis ###
avg_logFC <- 2
head(dropseq.combined@meta.data)
#dir.create('degs')

markers <- FindAllMarkers(object = dropseq.combined, min.pct = 0.25, logfc.threshold = 0.25, only.pos = T)
markers %>% group_by(cluster) %>% top_n(2, avg_logFC)
write.table(markers, 'degs/markers.res0.7.txt', sep = '\t', quote = F, col.names = T, row.names = F) 
top10 <- markers %>% group_by(cluster) %>% top_n(10, avg_logFC)
DoHeatmap(object = dropseq.combined, features = top10$gene, angle = 90, size = 1, raster = T, draw.lines = F)
#ggsave('degs/markers.res0.7.pdf', units = 'cm', width = 35, height = 25)

grep(rownames(dropseq.combined), value = T, pattern = 'muc2')
VlnPlot(object = dropseq.combined, pt.size = 0, features = c('cd8a', 'zap70'))
DotPlot(object = dropseq.combined, features = c('cdh5', 'sele', 'kdrl')) + coord_flip() # 16 vascular endothelium
DotPlot(object = dropseq.combined, features = c('tuba7l', 'dnai1.1', 'kif17')) + coord_flip() # 18 multiciliated cell
DotPlot(object = dropseq.combined, features = c('agr2', 'itln1')) + coord_flip() # 13 mucin cells
DotPlot(object = dropseq.combined, features = c('tnfrsf11a', 'spry2')) + coord_flip() # 18 nephron 14? 16?
DotPlot(object = dropseq.combined, features = c('clcnk', 'hoxb8a')) + coord_flip() # 14 distal tubule
DotPlot(object = dropseq.combined, features = c('slc26a1', 'slc5a12', 'slc13a3')) + coord_flip() # 17 proximal tubule
DotPlot(object = dropseq.combined, features = c('osr1', 'six2a', 'eya4')) + coord_flip() # 23 kidney progenitors
DotPlot(object = dropseq.combined, features = c('runx1t1', 'meis1b', 'tal1')) + coord_flip() # 21 HSCs
DotPlot(object = dropseq.combined, features = c('epor', 'sptb', 'ba1')) + coord_flip() # 2, 7, 12 Erythroid
DotPlot(object = dropseq.combined, features = c('mpx', 'lyz', 'cpa5')) + coord_flip() # 0, 4, 6, 11 Neutrophil
DotPlot(object = dropseq.combined, features = c('mpeg1.2', 'grn1', 'marco', 'grna', 'mpeg1.1', 'cts2.2')) + coord_flip() # 1, 5, 10, 20, 22 Macrophage
DotPlot(object = dropseq.combined, features = c('cd79a', 'pax5', 'igl3c3')) + coord_flip() # 15, 19 B cell
DotPlot(object = dropseq.combined, features = c('cd8a', 'zap70', 'cd4.1', 'jak3', 'il2rb', 'nkl.2', 'nkl.3', 'nkl.4', 'ccl33.3')) + coord_flip() #3, 8, 9 T/NK


#  0  1  2  3  4
#  5  6  7  8  9
# 10 11 12 13 14
# 15 16 17 18 19
# 20 21 22 23

celltypes <- c('Neutrophil', 'Macrophage', 'Erythroid', 'NK/T cell', 'Neutrophil',
               'Macrophage', 'Neutrophil', 'Erythroid', 'NK/T cell', 'NK/T cell',
               'Macrophage', 'Neutrophil', 'Erythroid', 'Mucin cells', 'Distal tubule',
               'B cell', 'Vascular endothelium', 'Proximal tubule', 'Nephron epithelium', 'B cell',
               'Macrophage', 'HSCs', 'Macrophage', 'Kidney progenitors')

dropseq.combined@meta.data$celltype <- mapvalues(dropseq.combined@meta.data$RNA_snn_res.0.7, from = levels(dropseq.combined@meta.data$RNA_snn_res.0.7), to = celltypes)
dropseq.combined@meta.data$celltype <- factor(dropseq.combined@meta.data$celltype, 
                                              levels = c("HSCs", "Erythroid", "Neutrophil", "Macrophage", "B cell", "NK/T cell", 
                                                         "Kidney progenitors", "Proximal tubule", "Distal tubule", "Nephron epithelium", "Mucin cells", "Vascular endothelium"))
Idents(dropseq.combined) <- 'celltype'
dropseq.combined@meta.data$seurat_clusters <- NULL
head(dropseq.combined@meta.data)
#saveRDS(dropseq.combined, 'tmp/zebrafish.harmony.Rds')

DimPlot(dropseq.combined, label = T, pt.size = 1, reduction = 'tsne')
ggsave('res/tsne.1_2.celltype.pdf', units = 'cm', width = 15.5, height = 10)
DimPlot(dropseq.combined, label = T, dims = c(1,3), pt.size = 1, reduction = 'tsne')
ggsave('res/tsne.1_3.celltype.pdf', units = 'cm', width = 15.5, height = 8.5)
DimPlot(dropseq.combined, label = T, dims = c(2,3), pt.size = 1, reduction = 'tsne')
ggsave('res/tsne.2_3.celltype.pdf', units = 'cm', width = 15.5, height = 8.5)

DimPlot(dropseq.combined, reduction = "umap", pt.size = 1, label = T)
ggsave('res/umap.1_2.min_dist_0.3.celltype.pdf', units = 'cm', width = 15.5, height = 10)
DimPlot(dropseq.combined, reduction = "umap", pt.size = 1, dims = c(1, 3), label = T)
ggsave('res/umap.1_3.min_dist_0.3.celltype.pdf', units = 'cm', width = 15.5, height = 10)
DimPlot(dropseq.combined, reduction = "umap", pt.size = 1, dims = c(2, 3), label = T)
ggsave('res/umap.2_3.min_dist_0.3.celltype.pdf', units = 'cm', width = 15.5, height = 10)


### Blood cells ###
bloodcells <- subset(dropseq.combined, celltype %in% c("HSCs", "Erythroid", "Neutrophil", "Macrophage", "B cell", "NK/T cell"))
bloodcells@meta.data <- droplevels(bloodcells@meta.data)
#saveRDS(bloodcells, 'tmp/zebrafish.harmony.bloodcells.Rds')

levels(Idents(bloodcells))
avg_logFC <- 2
markers <- FindAllMarkers(object = bloodcells, min.pct = 0.25, logfc.threshold = 0.25, only.pos = T, test.use = 'MAST')
markers %>% group_by(cluster) %>% top_n(2, avg_logFC)
write.table(markers, 'degs/markers_MAST.bloodcells.celltype.txt', sep = '\t', quote = F, col.names = T, row.names = F) 
top10 <- markers %>% group_by(cluster) %>% top_n(10, avg_logFC)
DoHeatmap(object = bloodcells, features = top10$gene, angle = 90, size = 5, raster = T, draw.lines = F)
ggsave('degs/markers_MAST.bloodcells.celltype.pdf', units = 'cm', width = 20, height = 15)

DimPlot(bloodcells, label = T, pt.size = 1, reduction = 'tsne')
ggsave('res/tsne.1_2.celltype.bloodcells.pdf', units = 'cm', width = 14, height = 10)
DimPlot(bloodcells, label = T, dims = c(1,3), pt.size = 1, reduction = 'tsne')
ggsave('res/tsne.1_3.celltype.bloodcells.pdf', units = 'cm', width = 14, height = 8.5)
DimPlot(bloodcells, label = T, dims = c(2,3), pt.size = 1, reduction = 'tsne')
ggsave('res/tsne.2_3.celltype.bloodcells.pdf', units = 'cm', width = 14, height = 8.5)

DimPlot(bloodcells, reduction = "umap", pt.size = 1, label = T)
ggsave('res/umap.1_2.min_dist_0.3.celltype.bloodcells.pdf', units = 'cm', width = 13.5, height = 10)
DimPlot(bloodcells, reduction = "umap", pt.size = 1, dims = c(1, 3), label = T)
ggsave('res/umap.1_3.min_dist_0.3.celltype.bloodcells.pdf', units = 'cm', width = 13.5, height = 10)
DimPlot(bloodcells, reduction = "umap", pt.size = 1, dims = c(2, 3), label = T)
ggsave('res/umap.2_3.min_dist_0.3.celltype.bloodcells.pdf', units = 'cm', width = 13.5, height = 10)



###
library(Seurat)
library(ggplot2)
library(plyr)
library(dplyr)

bloodcells <- readRDS('tmp/zebrafish.harmony.bloodcells.Rds')
Idents(bloodcells) <- 'celltype'
head(bloodcells@meta.data)
levels(bloodcells@meta.data$celltype)

fish_exprs <- as.matrix(GetAssayData(bloodcells, slot = 'data'))
fish_exprs <- exp(fish_exprs) -1; head(colSums(fish_exprs))
fish_exprs[1:4, 1:4]


### Pseudobulk ###
genes <- rownames(bloodcells); length(genes)
celltype_pseudo <- data.frame(matrix(nrow = length(genes), ncol = length(levels(bloodcells@meta.data$celltype))))
rownames(celltype_pseudo) <- genes
colnames(celltype_pseudo) <- levels(bloodcells@meta.data$celltype)
head(celltype_pseudo)

for (ct in levels(bloodcells@meta.data$celltype)) {
  ct_cells <- rownames(subset(bloodcells@meta.data, celltype == ct))
  
  ct_pseudo <- fish_exprs[, ct_cells]
  ct_pseudo <- rowMeans(ct_pseudo)
  celltype_pseudo[, ct] <- ct_pseudo
}

celltype_pseudo$sums <- rowSums(celltype_pseudo)
celltype_pseudo <- subset(celltype_pseudo, sums != 0)
celltype_pseudo$sums <- NULL
head(celltype_pseudo); dim(celltype_pseudo)
#saveRDS(celltype_pseudo, 'tmp/celltype.pseudo.Rds')


### Pseudo-cell (10 cells) transformation ###
pseudocell_total <- data.frame(matrix(ncol = 0, nrow = length(rownames(bloodcells))))
rownames(pseudocell_total) <- rownames(bloodcells)

for (ct in levels(bloodcells@meta.data$celltype)) {
  tmpObj <- subset(bloodcells, celltype == ct)
  tmpObj_exprs <- as.matrix(GetAssayData(tmpObj, slot = 'data', assay = 'RNA'))
  
  bcs <- rownames(tmpObj@meta.data)
  rounds <- trunc(length(bcs)/10)
  
  ct_pseudocells <- data.frame(matrix(ncol = rounds, nrow = length(rownames(bloodcells))))
  colnames(ct_pseudocells) <- unlist(lapply(c(1:rounds), function (x) paste0(c(as.character(ct), x), collapse = ' pc') ))
  rownames(ct_pseudocells) <- rownames(bloodcells)
  
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
saveRDS(pseudocell_total, 'tmp/celltype.pseudocell10.Rds')


#
library(Rtsne)

bloodcells <- FindVariableFeatures(object = bloodcells, selection.method = "vst", nfeatures = 2000)
#saveRDS(bloodcells, 'tmp/zebrafish.harmony.bloodcells.Rds')
varigene_mat <- t(pseudocell_total[VariableFeatures(bloodcells), ])

tsne <- Rtsne(as.matrix(varigene_mat), 
              dims = 3, check_duplicates = F, theta = 0.5, pca_scale = T)

scores <- as.data.frame(tsne$Y)
head(scores); nrow(scores)
rownames(scores) <- colnames(pseudocell_total)
colnames(scores) <- c('tSNE1', 'tSNE2', 'tSNE3')
scores$Celltype <- unlist( lapply(rownames(scores), function (x) unlist(strsplit(as.character(x), split = ' pc'))[1] ) )
scores$Celltype <- factor(scores$Celltype, levels = levels(bloodcells@meta.data$celltype))
head(scores)

ggplot(scores, aes(x = tSNE1, y = tSNE2, colour = Celltype)) + 
  geom_point(size = 2) + 
  scale_colour_manual(name="", values = c("HSCs"="grey60", "Erythroid"="#005DBC", "Neutrophil"="#b65aa9", 
                                          "Macrophage"="#feb300", "B cell"="#006400", "NK/T cell"="#EE0000")) + 
  theme_bw() + 
  theme(axis.text = element_text(colour = 'black'),
        axis.ticks = element_line(colour = 'black'),
        panel.grid = element_blank()) +
  labs(title =  '')
#ggsave('tsne.pseudocell10.pdf', units = 'cm', width = 11, height = 8)


exprs <- exp(as.matrix(GetAssayData(bloodcells, slot = 'data'))) -1
head(colSums(exprs))
exprs <- exprs[VariableFeatures(bloodcells), ]
tsne <- Rtsne(as.matrix(t(exprs)), 
              dims = 3, check_duplicates = F, theta = 0.5, pca_scale = T)

scores <- as.data.frame(tsne$Y)
head(scores); nrow(scores)
rownames(scores) <- colnames(exprs)
colnames(scores) <- c('tSNE1', 'tSNE2', 'tSNE3')
scores$Celltype <- bloodcells@meta.data$celltype
scores$Celltype <- factor(scores$Celltype, levels = levels(bloodcells@meta.data$celltype))
head(scores)

ggplot(scores, aes(x = tSNE1, y = tSNE2, colour = Celltype)) + 
  geom_point(size = 1) + 
  scale_colour_manual(name="", values = c("HSCs"="grey60", "Erythroid"="#005DBC", "Neutrophil"="#b65aa9", 
                                          "Macrophage"="#feb300", "B cell"="#006400", "NK/T cell"="#EE0000")) + 
  theme_bw() + 
  theme(axis.text = element_text(colour = 'black'),
        axis.ticks = element_line(colour = 'black'),
        panel.grid = element_blank()) +
  labs(title =  '')
#ggsave('tsne.allcells.pdf', units = 'cm', width = 11, height = 8)



###
library(Seurat)
library(ggplot2)
library(plyr)
library(dplyr)

bloodcells <- readRDS('tmp/zebrafish.harmony.bloodcells.Rds')
head(bloodcells@meta.data)
levels(bloodcells@meta.data$celltype)

dir.create('degs')
degmarkers <- FindAllMarkers(object = bloodcells, min.pct = 0.25, logfc.threshold = 0.25, only.pos = T, test.use = 'MAST')
write.table(degmarkers, 'degs/celltype.txt', sep = '\t', quote = F, col.names = T, row.names = F) 
saveRDS(degmarkers, 'degs/celltype.Rds')

top10 <- degmarkers %>% group_by(cluster) %>% top_n(10, avg_logFC)
dehm <- DoHeatmap(object = bloodcells, features = top10$gene, angle = 90, size = 3, raster = T, draw.lines = F)
ggsave('degs/celltype.pdf', units = 'cm', width = 30, height = 20)



###
library(Seurat)
library(plyr)
library(dplyr)
library(ggplot2)
library(RColorBrewer)

zebrafish <- readRDS('tmp/zebrafish.harmony.bloodcells.Rds')
head(zebrafish@meta.data)

zebrafish <- JackStraw(object = zebrafish, num.replicate = 100, dims = 40)
zebrafish <- ScoreJackStraw(object = zebrafish, dims = 1:40)
JackStrawPlot(zebrafish, dims = 1:40) # 22 PCs 
ggsave('stats/2-4.pca.JackStrawPlot.mitoCut.bloodcells.pdf', units = 'cm', width = 20, height = 12)

zebrafish <- RunUMAP(zebrafish, dims=1:22, reduction = "harmony", reduction.key='UMAP', n.components=3, min.dist=0.3, seed.use = 21060702)

DimPlot(zebrafish, pt.size = 1, reduction = 'umap') + 
  scale_color_manual(values = colorRampPalette(brewer.pal(9, "Spectral"))(length(levels(zebrafish@meta.data$celltype))) ) + 
  theme_void() + 
  theme(legend.position = 'none')
ggsave('res/new.umap.seed21060702.1_2.celltype_refined.void.png', units = 'cm', width = 8, height = 8)
DimPlot(zebrafish, pt.size = 1, reduction = 'umap', label = T, label.size = 2) + 
  scale_color_manual(values = colorRampPalette(brewer.pal(9, "Spectral"))(length(levels(zebrafish@meta.data$celltype))) )
ggsave('res/new.umap.seed21060702.1_2.celltype_refined.pdf', units = 'cm', width = 11.5, height = 8)

saveRDS(zebrafish, 'tmp/zebrafish.harmony.bloodcells.Rds')


DotPlot(zebrafish, features = c('mmp14b', 'slc30a1a', 'laptm4a', 'skp1', 'npepl1', 'hspa5')) +
  labs(x = '', y = '') + coord_flip() + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5))
ggsave('degs/celltype_unique.pm120.pdf', units = 'cm', width = 12, height = 10)



