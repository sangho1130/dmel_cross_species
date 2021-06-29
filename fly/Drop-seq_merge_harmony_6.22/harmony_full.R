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


exprs <- '/home/sangho/2020_newbuild/projects/drosophila/projectX_metabolism/table/merged/merged.count.6.22.Rds'
labels <- '/home/sangho/2020_newbuild/projects/drosophila/projectX_metabolism/table/merged/merged.label.Rds'

data <- readRDS(exprs)
labelData <- readRDS(labels)
colnames(labelData)[1] <- 'dataset'
head(labelData)

dropseq <- CreateSeuratObject(counts = data, project = "Dmel_harmony")
dropseq <- AddMetaData(object = dropseq, metadata = labelData)
nrow(dropseq@meta.data) #  cells

plt <- ggplot(dropseq@meta.data, aes(Library, nFeature_RNA)) + geom_jitter(size = 0.25) + 
  theme_bw() + 
  theme(text = element_text(size = 7), axis.text = element_text(size = 7), panel.grid = element_blank(), axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5));plt
ggsave('stats/2-2.stats.nFeature_RNA.byLibrary.pdf', units = 'cm', width = 30, height = 8)
plt <- ggplot(dropseq@meta.data, aes(Library, nCount_RNA)) + geom_jitter(size = 0.25) + 
  theme_bw() + 
  theme(text = element_text(size = 7), axis.text = element_text(size = 7), panel.grid = element_blank(), axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5));plt
ggsave('stats/2-2.stats.nCount_RNA.byLibrary.pdf', units = 'cm', width = 30, height = 8)

plt <- ggplot(dropseq@meta.data, aes(Library, percent.mt)) + geom_jitter(size = 0.25) + 
  theme_bw() +
  theme(text = element_text(size = 7), axis.text = element_text(size = 7), panel.grid = element_blank(), axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5));plt
ggsave('stats/2-2.stats.percent.mt.byLibrary.pdf', units = 'cm', width = 30, height = 8)

plt <- ggplot(dropseq@meta.data, aes(percent.mt, col = Library)) + geom_density() +
  theme_bw() + theme(text = element_text(size = 7), axis.text = element_text(size = 7), panel.grid = element_blank());plt
ggsave('stats/2-2.stats.percent.mt.density_byLibrary.pdf', units = 'cm', width = 24, height = 12)
plt <- ggplot(dropseq@meta.data, aes(nCount_RNA, col = Library)) + geom_density() +
  theme_bw() + theme(text = element_text(size = 7), axis.text = element_text(size = 7), panel.grid = element_blank());plt
ggsave('stats/2-2.stats.nCount_RNA.density_byLibrary.pdf', units = 'cm', width = 12, height = 6)
plt <- ggplot(dropseq@meta.data, aes(nFeature_RNA, col = Library)) + geom_density() +
  theme_bw() + theme(text = element_text(size = 7), axis.text = element_text(size = 7), panel.grid = element_blank());plt
ggsave('stats/2-2.stats.nFeature_RNA.density_byLibrary.pdf', units = 'cm', width = 12, height = 6)

plot1 <- FeatureScatter(dropseq, feature1 = "nFeature_RNA", feature2 = "percent.mt", group.by = 'Library')
plot2 <- FeatureScatter(dropseq, feature1 = "nCount_RNA", feature2 = "percent.mt", group.by = 'Library')
CombinePlots(plots = list(plot1, plot2))
ggsave('stats/2-2.stats.percent.mt.compare.pdf', units = 'cm', width = 40, height = 10)


### The standard workflow ###
dropseq.combined <- NormalizeData(dropseq, normalization.method = "LogNormalize", scale.factor = 10000)
dropseq.combined <- FindVariableFeatures(object = dropseq.combined, selection.method = "vst", nfeatures = 2000)

dropseq.combined <- ScaleData(object = dropseq.combined, vars.to.regress = c('Library', 'nCount_RNA'))
dropseq.combined <- RunPCA(object = dropseq.combined, npcs = 80)
dropseq.combined <- JackStraw(object = dropseq.combined, num.replicate = 100, dims = 80)
dropseq.combined <- ScoreJackStraw(object = dropseq.combined, dims = 1:80)
JackStrawPlot(dropseq.combined, dims = 1:80) # 51 PCs 
ggsave('stats/2-4.pca.JackStrawPlot.mitoCut.pdf', units = 'cm', width = 30, height = 12)
saveRDS(dropseq.combined, 'tmp/dropseq.combined_2.Rds')

dropseq.combined <- RunHarmony(dropseq.combined, group.by.vars = "dataset")
saveRDS(dropseq.combined, 'tmp/dropseq.combined.harmony.Rds')

###
dropseq.combined <- readRDS('tmp/dropseq.combined.harmony.Rds')
genes <- c("Antp", "Dl", "NimB3", "IM18", "Ance", "Hml", "NimC1", "mthl4", "PPO1")
### UMAP ###
#for (seed in c(210119250:210119259)){
for (seed in c(210120150:210120199)){ # second round
  dropseq.combined <- RunUMAP(dropseq.combined, dims=1:51, reduction = "harmony", reduction.key='UMAP', n.components=3, min.dist=0.2, seed.use = seed)
  FeaturePlot(dropseq.combined, features = genes, cols = c("grey","red"), reduction = "umap")
  ggsave(paste0(c('umap/compare/umap.1_2.markers.seed', seed, '.png'), collapse = ''), units = 'cm', width = 35, height = 30)
} # 210120153 (min.dist 0.2)

### tSNE ###
#for (seed in c(210119200:210119209)){
for (seed in c(210120100:210120149)){ # second round
  dropseq.combined <- RunTSNE(dropseq.combined, dims=1:51, reduction = "harmony", reduction.key='tSNE', dim.embed=3, seed.use = seed)
  FeaturePlot(dropseq.combined, features = genes, cols = c("grey","red"), reduction = "tsne")
  ggsave(paste0(c('tsne/compare/tsne.1_2.markers.seed', seed, '.png'), collapse = ''), units = 'cm', width = 35, height = 30)
} # 210119203

### UMAP fix ###
dropseq.combined <- RunUMAP(dropseq.combined, dims=1:51, reduction = "harmony", reduction.key='UMAP', n.components=3, min.dist=0.2, seed.use = 210120153)
FeaturePlot(dropseq.combined, features = genes, cols = c("grey","red"), reduction = "umap")
ggsave('umap/umap.combined.1_2.markers.seed210120153.png', units = 'cm', width = 35, height = 30)

DimPlot(object=dropseq.combined, group.by = "dataset",  reduction = 'umap')
ggsave('umap/umap.combined.1_2.dataset.seed210120153.pdf', units = 'cm', width = 15, height = 10)
DimPlot(object=dropseq.combined, group.by = "Library", reduction = 'umap')
ggsave('umap/umap.combined.1_2.libraries.seed210120153.pdf', units = 'cm', width = 17, height = 10)
FeaturePlot(dropseq.combined, features = 'nFeature_RNA', cols = c("grey","red"), reduction = "umap")
ggsave('umap/umap.combined.1_2.nFeature_RNA.seed210120153.pdf', units = 'cm', width = 12, height = 10)
FeaturePlot(dropseq.combined, features = 'nCount_RNA', cols = c("grey","red"), reduction = "umap")
ggsave('umap/umap.combined.1_2.nCount_RNA.seed210120153.pdf', units = 'cm', width = 12, height = 10)
FeaturePlot(dropseq.combined, features = 'percent.mt', cols = c("grey","red"), reduction = "umap")
ggsave('umap/umap.combined.1_2.percent.mt.seed210120153.pdf', units = 'cm', width = 11, height = 10)
DimPlot(object=dropseq.combined, group.by = "origin",  reduction = 'umap')
ggsave('umap/umap.combined.1_2.origin.seed210120153.pdf', units = 'cm', width = 14, height = 10)
DimPlot(object=dropseq.combined, group.by = "condition",  reduction = 'umap')
ggsave('umap/umap.combined.1_2.condition.seed210120153.pdf', units = 'cm', width = 14, height = 10)

### tSNE fix ###
dropseq.combined <- RunTSNE(dropseq.combined, dims=1:51, reduction = "harmony", reduction.key='tSNE', dim.embed=3, seed.use = 210119203)
FeaturePlot(dropseq.combined, features = genes, cols = c("grey","red"), reduction = "tsne")
ggsave('tsne/tsne.combined.1_2.markers.seed210119203.png', units = 'cm', width = 35, height = 30)

DimPlot(object=dropseq.combined, group.by = "dataset", reduction = 'tsne')
ggsave('tsne/tsne.combined.1_2.dataset.seed210119203.pdf', units = 'cm', width = 15, height = 10)
DimPlot(object=dropseq.combined, group.by = "Library", reduction = 'tsne')
ggsave('tsne/tsne.combined.1_2.libraries.seed210119203.pdf', units = 'cm', width = 17, height = 10)
FeaturePlot(dropseq.combined, features = 'nFeature_RNA', cols = c("grey","red"), reduction = "tsne")
ggsave('tsne/tsne.combined.1_2.nFeature_RNA.seed210119203.pdf', units = 'cm', width = 12, height = 10)
FeaturePlot(dropseq.combined, features = 'nCount_RNA', cols = c("grey","red"), reduction = "tsne")
ggsave('tsne/tsne.combined.1_2.nCount_RNA.seed210119203.pdf', units = 'cm', width = 12, height = 10)
FeaturePlot(dropseq.combined, features = 'percent.mt', cols = c("grey","red"), reduction = "tsne")
ggsave('tsne/tsne.combined.1_2.percent.mt.seed210119203.pdf', units = 'cm', width = 11, height = 10)
DimPlot(object=dropseq.combined, group.by = "origin",  reduction = 'tsne')
ggsave('tsne/tsne.combined.1_2.origin.seed210119203.pdf', units = 'cm', width = 14, height = 10)
DimPlot(object=dropseq.combined, group.by = "condition",  reduction = 'tsne')
ggsave('tsne/tsne.combined.1_2.condition.seed210119203.pdf', units = 'cm', width = 14, height = 10)

saveRDS(dropseq.combined, 'tmp/dropseq.combined.harmony.Rds')


### Filter non-hematopoietic cells 
library(Seurat)
library(ggplot2)
library(plyr)
library(dplyr)
library(harmony)
dir.create('stats_flt')
dir.create('umap/compare_flt/')
dir.create('tsne/compare_flt/')

dropseq.combined <- readRDS('tmp/dropseq.combined.harmony.Rds')
dropseq.combined <- subset(dropseq.combined, cells = rownames(subset(dropseq.combined@meta.data, !celltype %in% c('DV', 'RG', 'Neurons')))) # additional 42 PSC from circulationg PI24 were removed
rmcells <- rownames(subset(dropseq.combined@meta.data, origin == 'Circulation' & celltype == 'PSC'))
dropseq.combined <- subset(dropseq.combined, cells = setdiff(rownames(dropseq.combined@meta.data), rmcells))

dropseq.combined@meta.data <- droplevels(dropseq.combined@meta.data)
head(dropseq.combined@meta.data)

dropseq.combined <- NormalizeData(dropseq.combined, normalization.method = "LogNormalize", scale.factor = 10000)
dropseq.combined <- FindVariableFeatures(object = dropseq.combined, selection.method = "vst", nfeatures = 2000)
dropseq.combined <- ScaleData(object = dropseq.combined, vars.to.regress = c('Library', 'nCount_RNA'))
dropseq.combined <- RunPCA(object = dropseq.combined, npcs = 80)
dropseq.combined <- JackStraw(object = dropseq.combined, num.replicate = 100, dims = 80)
dropseq.combined <- ScoreJackStraw(object = dropseq.combined, dims = 1:80)
JackStrawPlot(dropseq.combined, dims = 1:80) # 50 PCs 
ggsave('stats_flt/2-4.pca.JackStrawPlot.mitoCut.pdf', units = 'cm', width = 30, height = 12)

dropseq.combined <- RunHarmony(dropseq.combined, group.by.vars = "dataset")
saveRDS(dropseq.combined, 'tmp/dropseq.combined.harmony_flt.Rds')

for (seed in c(210205150:210205199)){ 
  dropseq.combined <- RunUMAP(dropseq.combined, dims=1:50, reduction = "harmony", reduction.key='UMAP', n.components=3, min.dist=0.25, spread = 4, seed.use = seed)
  DimPlot(object=dropseq.combined, group.by = "celltype",  reduction = 'umap') +
    scale_color_manual(values = c('#F15FA5','#1F7EB2','#A70D0C','#F0A041','#25A8B0','#a4a4a4','#1a1a1a'))
  ggsave(paste0(c('umap/compare_flt/umap.1_2.celltype.seed', seed, '.png'), collapse = ''), units = 'cm', width = 14, height = 10)
  DimPlot(object=dropseq.combined, group.by = "celltype",  reduction = 'umap', dims = c(1,3)) +
    scale_color_manual(values = c('#F15FA5','#1F7EB2','#A70D0C','#F0A041','#25A8B0','#a4a4a4','#1a1a1a'))
  ggsave(paste0(c('umap/compare_flt/umap.1_3.celltype.seed', seed, '.png'), collapse = ''), units = 'cm', width = 14, height = 10)
  DimPlot(object=dropseq.combined, group.by = "celltype",  reduction = 'umap', dims = c(2,3)) +
    scale_color_manual(values = c('#F15FA5','#1F7EB2','#A70D0C','#F0A041','#25A8B0','#a4a4a4','#1a1a1a'))
  ggsave(paste0(c('umap/compare_flt/umap.2_3.celltype.seed', seed, '.png'), collapse = ''), units = 'cm', width = 14, height = 10)
} # 210202191 (0.2, spread=1), 210202353 (0.2, spread=2), 210202388 (0.2, spread=1.5), 210203360 (0.2, spread=2), 210205195 (0.25, spread=4)

for (seed in c(210202100:210202149)){
  dropseq.combined <- RunTSNE(dropseq.combined, dims=1:50, reduction = "harmony", reduction.key='tSNE', dim.embed=3, seed.use = seed)
  DimPlot(object=dropseq.combined, group.by = "celltype",  reduction = 'tsne') +
    scale_color_manual(values = c('#F15FA5','#1F7EB2','#A70D0C','#F0A041','#25A8B0','#a4a4a4','#1a1a1a'))
  ggsave(paste0(c('tsne/compare_flt/tsne.1_2.markers.seed', seed, '.png'), collapse = ''), units = 'cm', width = 14, height = 10)
} # 210202127

dropseq.combined <- RunUMAP(dropseq.combined, dims=1:50, reduction = "harmony", reduction.key='UMAP', n.components=3, min.dist=0.25, spread = 4, seed.use = 210205195)
dropseq.combined <- RunTSNE(dropseq.combined, dims=1:50, reduction = "harmony", reduction.key='tSNE', dim.embed=3, seed.use = 210202127)
saveRDS(dropseq.combined, 'tmp/dropseq.combined.harmony_flt.Rds')


DimPlot(object=dropseq.combined, group.by = "dataset",  reduction = 'umap')
ggsave('umap/umap.combined.1_2.dataset.seed210205195.pdf', units = 'cm', width = 15, height = 10)
DimPlot(object=dropseq.combined, group.by = "Library", reduction = 'umap')
ggsave('umap/umap.combined.1_2.libraries.seed210205195.pdf', units = 'cm', width = 17, height = 10)
FeaturePlot(dropseq.combined, features = 'nFeature_RNA', cols = c("grey","red"), reduction = "umap")
ggsave('umap/umap.combined.1_2.nFeature_RNA.seed210205195.pdf', units = 'cm', width = 12, height = 10)
FeaturePlot(dropseq.combined, features = 'nCount_RNA', cols = c("grey","red"), reduction = "umap")
ggsave('umap/umap.combined.1_2.nCount_RNA.seed210205195.pdf', units = 'cm', width = 12, height = 10)
FeaturePlot(dropseq.combined, features = 'percent.mt', cols = c("grey","red"), reduction = "umap")
ggsave('umap/umap.combined.1_2.percent.mt.seed210205195.pdf', units = 'cm', width = 11, height = 10)
DimPlot(object=dropseq.combined, group.by = "origin",  reduction = 'umap')
ggsave('umap/umap.combined.1_2.origin.seed210205195.pdf', units = 'cm', width = 14, height = 10)
DimPlot(object=dropseq.combined, group.by = "condition",  reduction = 'umap')
ggsave('umap/umap.combined.1_2.condition.seed210205195.pdf', units = 'cm', width = 14, height = 10)

DimPlot(object=dropseq.combined, group.by = "dataset",  reduction = 'umap', dims = c(1,3))
ggsave('umap/umap.combined.1_3.dataset.seed210205195.pdf', units = 'cm', width = 15, height = 10)
DimPlot(object=dropseq.combined, group.by = "Library", reduction = 'umap', dims = c(1,3))
ggsave('umap/umap.combined.1_3.libraries.seed210205195.pdf', units = 'cm', width = 17, height = 10)
FeaturePlot(dropseq.combined, features = 'nFeature_RNA', cols = c("grey","red"), reduction = "umap", dims = c(1,3))
ggsave('umap/umap.combined.1_3.nFeature_RNA.seed210205195.pdf', units = 'cm', width = 12, height = 10)
FeaturePlot(dropseq.combined, features = 'nCount_RNA', cols = c("grey","red"), reduction = "umap", dims = c(1,3))
ggsave('umap/umap.combined.1_3.nCount_RNA.seed210205195.pdf', units = 'cm', width = 12, height = 10)
FeaturePlot(dropseq.combined, features = 'percent.mt', cols = c("grey","red"), reduction = "umap", dims = c(1,3))
ggsave('umap/umap.combined.1_3.percent.mt.seed210205195.pdf', units = 'cm', width = 11, height = 10)
DimPlot(object=dropseq.combined, group.by = "origin",  reduction = 'umap', dims = c(1,3))
ggsave('umap/umap.combined.1_3.origin.seed210205195.pdf', units = 'cm', width = 14, height = 10)
DimPlot(object=dropseq.combined, group.by = "condition",  reduction = 'umap', dims = c(1,3))
ggsave('umap/umap.combined.1_3.condition.seed210205195.pdf', units = 'cm', width = 14, height = 10)



DimPlot(object=dropseq.combined, group.by = "dataset", reduction = 'tsne')
ggsave('tsne/tsne.combined.1_2.dataset.seed210202127.pdf', units = 'cm', width = 15, height = 10)
DimPlot(object=dropseq.combined, group.by = "Library", reduction = 'tsne')
ggsave('tsne/tsne.combined.1_2.libraries.seed210202127.pdf', units = 'cm', width = 17, height = 10)
FeaturePlot(dropseq.combined, features = 'nFeature_RNA', cols = c("grey","red"), reduction = "tsne")
ggsave('tsne/tsne.combined.1_2.nFeature_RNA.seed210202127.pdf', units = 'cm', width = 12, height = 10)
FeaturePlot(dropseq.combined, features = 'nCount_RNA', cols = c("grey","red"), reduction = "tsne")
ggsave('tsne/tsne.combined.1_2.nCount_RNA.seed210202127.pdf', units = 'cm', width = 12, height = 10)
FeaturePlot(dropseq.combined, features = 'percent.mt', cols = c("grey","red"), reduction = "tsne")
ggsave('tsne/tsne.combined.1_2.percent.mt.seed210202127.pdf', units = 'cm', width = 11, height = 10)
DimPlot(object=dropseq.combined, group.by = "origin",  reduction = 'tsne')
ggsave('tsne/tsne.combined.1_2.origin.seed210202127.pdf', units = 'cm', width = 14, height = 10)
DimPlot(object=dropseq.combined, group.by = "condition",  reduction = 'tsne')
ggsave('tsne/tsne.combined.1_2.condition.seed210202127.pdf', units = 'cm', width = 14, height = 10)

DimPlot(object=dropseq.combined, group.by = "dataset", reduction = 'tsne', dims = c(1,3))
ggsave('tsne/tsne.combined.1_3.dataset.seed210202127.pdf', units = 'cm', width = 15, height = 10)
DimPlot(object=dropseq.combined, group.by = "Library", reduction = 'tsne', dims = c(1,3))
ggsave('tsne/tsne.combined.1_3.libraries.seed210202127.pdf', units = 'cm', width = 17, height = 10)
FeaturePlot(dropseq.combined, features = 'nFeature_RNA', cols = c("grey","red"), reduction = "tsne", dims = c(1,3))
ggsave('tsne/tsne.combined.1_3.nFeature_RNA.seed210202127.pdf', units = 'cm', width = 12, height = 10)
FeaturePlot(dropseq.combined, features = 'nCount_RNA', cols = c("grey","red"), reduction = "tsne", dims = c(1,3))
ggsave('tsne/tsne.combined.1_3.nCount_RNA.seed210202127.pdf', units = 'cm', width = 12, height = 10)
FeaturePlot(dropseq.combined, features = 'percent.mt', cols = c("grey","red"), reduction = "tsne", dims = c(1,3))
ggsave('tsne/tsne.combined.1_3.percent.mt.seed210202127.pdf', units = 'cm', width = 11, height = 10)
DimPlot(object=dropseq.combined, group.by = "origin",  reduction = 'tsne', dims = c(1,3))
ggsave('tsne/tsne.combined.1_3.origin.seed210202127.pdf', units = 'cm', width = 14, height = 10)
DimPlot(object=dropseq.combined, group.by = "condition",  reduction = 'tsne', dims = c(1,3))
ggsave('tsne/tsne.combined.1_3.condition.seed210202127.pdf', units = 'cm', width = 14, height = 10)
