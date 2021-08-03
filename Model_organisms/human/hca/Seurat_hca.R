library(Seurat)
library(ggplot2)
library(plyr)
library(dplyr)

dir.create('tmp')
dir.create('stats')
dir.create('umap')
dir.create('umap/compare')
dir.create('degs')

counts <- readRDS('allcells.counts.BM.Rds')
label <- readRDS('label.BM.Rds')
label$barcode <- NULL
colnames(label)[1] <- 'Library'
head(label)
identical(rownames(label), colnames(counts))


### The standard workflow ###
human_hca <- CreateSeuratObject(counts = counts, project = "HCA")
human_hca <- AddMetaData(object = human_hca, metadata = label)
head(human_hca@meta.data); nrow(human_hca@meta.data) # 290861 cells
human_hca@meta.data$Library <- factor(human_hca@meta.data$Library, levels = c('donor1','donor2','donor3','donor4','donor5','donor6','donor7','donor8'))
summary(human_hca@meta.data)

plt <- ggplot(human_hca@meta.data, aes(Library, nFeature_RNA)) + geom_jitter(size = 0.25) + 
  theme_bw() + 
  theme(text = element_text(size = 7), axis.text = element_text(size = 7), panel.grid = element_blank(), axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5));plt
ggsave('stats/2-2.stats.human_hca.nFeature_RNA.byLibrary.png', units = 'cm', width = 6, height = 6)
plt <- ggplot(human_hca@meta.data, aes(Library, nCount_RNA)) + geom_jitter(size = 0.25) + 
  theme_bw() + geom_hline(yintercept = 80000, col = 'red2') +
  theme(text = element_text(size = 7), axis.text = element_text(size = 7), panel.grid = element_blank(), axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5));plt
ggsave('stats/2-2.stats.human_hca.nCount_RNA.byLibrary.png', units = 'cm', width = 6, height = 6)

human_hca <- subset(human_hca, nCount_RNA < 80000)
human_hca <- subset(human_hca, nFeature_RNA >=  500)
plt <- ggplot(human_hca@meta.data, aes(Library, nFeature_RNA)) + geom_jitter(size = 0.25) + 
  theme_bw() + 
  theme(text = element_text(size = 7), axis.text = element_text(size = 7), panel.grid = element_blank(), axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5));plt
ggsave('stats/2-3.stats.human_hca.nFeature_RNA.byLibrary.png', units = 'cm', width = 6, height = 6)
plt <- ggplot(human_hca@meta.data, aes(Library, nCount_RNA)) + geom_jitter(size = 0.25) + 
  theme_bw() + 
  theme(text = element_text(size = 7), axis.text = element_text(size = 7), panel.grid = element_blank(), axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5));plt
ggsave('stats/2-3.stats.human_hca.nCount_RNA.byLibrary.png', units = 'cm', width = 6, height = 6)

plt <- ggplot(human_hca@meta.data, aes(Library, nCount_RNA)) + geom_jitter(size = 0.25) + 
  theme_bw() + ylim(0, 5000) +  
  theme(text = element_text(size = 7), axis.text = element_text(size = 7), panel.grid = element_blank(), axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5));plt
ggsave('stats/2-3.stats.human_hca.nCount_RNA.byLibrary.lower.png', units = 'cm', width = 6, height = 6)

saveRDS(human_hca, 'tmp/human_hca1.Rds') # 260764 cells


### Count & MT threshold ###
head(human_hca@meta.data)
split_human_hca <- SplitObject(human_hca, split.by='Library')
for (i in c(1:8)) {
  umithre <- as.integer(mean(split_human_hca[[i]]@meta.data$nCount_RNA) + 2*sd(split_human_hca[[i]]@meta.data$nCount_RNA))
  split_human_hca[[i]] <- subset(split_human_hca[[i]], subset = nCount_RNA < umithre)
  
  split_human_hca[[i]][["percent.mt"]] <- PercentageFeatureSet(object = split_human_hca[[i]], pattern = "^MT-")
  split_human_hca[[i]] <- subset(split_human_hca[[i]], subset = percent.mt <= 10)
  
  split_human_hca[[i]] <- NormalizeData(split_human_hca[[i]], normalization.method = "LogNormalize", scale.factor = 10000)
  split_human_hca[[i]] <- FindVariableFeatures(object = split_human_hca[[i]], selection.method = "vst", nfeatures = 2000)
}


### Alignment ###
multi.anchors <- FindIntegrationAnchors(object.list = split_human_hca, dims = 1:30)
saveRDS(multi.anchors, 'tmp/multi.anchors.Rds')
human_hca <- IntegrateData(anchorset = multi.anchors, dims = 1:30)
head(human_hca@meta.data); nrow(human_hca@meta.data)
DefaultAssay(object = human_hca) <- "integrated"
saveRDS(human_hca, 'tmp/human_hca2.Rds')

plt <- ggplot(human_hca@meta.data, aes(Library, nFeature_RNA)) + geom_jitter(size = 0.25) + 
  theme_bw() + 
  theme(text = element_text(size = 7), axis.text = element_text(size = 7), panel.grid = element_blank(), axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5));plt
ggsave('stats/2-4.stats.human_hca.nFeature_RNA.byLibrary.png', units = 'cm', width = 6, height = 6)
plt <- ggplot(human_hca@meta.data, aes(Library, nCount_RNA)) + geom_jitter(size = 0.25) + 
  theme_bw() + 
  theme(text = element_text(size = 7), axis.text = element_text(size = 7), panel.grid = element_blank(), axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5));plt
ggsave('stats/2-4.stats.human_hca.nCount_RNA.byLibrary.png', units = 'cm', width = 6, height = 6)


### Scale & PCA ###
#DefaultAssay(object = human_hca) <- "RNA"
human_hca <- readRDS('tmp/human_hca2.Rds')
human_hca <- FindVariableFeatures(object = human_hca, selection.method = "vst", nfeatures = 2000)

human_hca <- ScaleData(object = human_hca, vars.to.regress = c('Library', 'nCount_RNA'))
human_hca <- RunPCA(object = human_hca, npcs = 80)
human_hca <- JackStraw(object = human_hca, num.replicate = 100, dims = 80)
human_hca <- ScoreJackStraw(object = human_hca, dims = 1:80)
JackStrawPlot(human_hca, dims = 1:80) # 59 PCs 
ggsave('stats/integrated.2-4.pca.JackStrawPlot.mitoCut.pdf', units = 'cm', width = 25, height = 12)
saveRDS(human_hca, 'tmp/integrated.human_hca.alignment.Rds')


### UMAP ###
genes <- c("CD34", "CD3E", "NKG7", "CD79A", "TPSAB1", "LYZ", "C1QB", "FCGR3A", "COL1A1")
for (seed in c(210225150:210225179)){
  human_hca <- RunUMAP(human_hca, dims=1:59, reduction = "pca", reduction.key='UMAP', n.components=3, min.dist=0.3, seed.use = seed)
  FeaturePlot(human_hca, features = genes, cols = c("grey", "red"), pt.size = 0.5, reduction = "umap")
  ggsave(paste0(c('umap/compare/integrated.umap.1_2.markers.seed', seed, '.png'), collapse = ''), units = 'cm', width = 35, height = 30)
} # 210225166


### UMAP fix ###
head(human_hca@meta.data)
human_hca <- RunUMAP(human_hca, dims=1:59, reduction = "pca", reduction.key='UMAP', n.components=3, min.dist=0.3, seed.use = 210225166)
FeaturePlot(human_hca, features = genes, cols = c("grey", "red"), pt.size = 0.5, reduction = "umap")
ggsave('umap/integrated.umap.combined.1_2.markers.seed210225166.png', units = 'cm', width = 35, height = 30)
DimPlot(object=human_hca, group.by = "Library", reduction = 'umap')
ggsave('umap/integrated.umap.combined.1_2.libraries.seed210225166.png', units = 'cm', width = 13, height = 10)
FeaturePlot(human_hca, features = 'nFeature_RNA', cols = c("grey","red"), reduction = "umap")
ggsave('umap/integrated.umap.combined.1_2.nFeature_RNA.seed210225166.png', units = 'cm', width = 12, height = 10)
FeaturePlot(human_hca, features = 'nCount_RNA', cols = c("grey","red"), reduction = "umap")
ggsave('umap/integrated.umap.combined.1_2.nCount_RNA.seed210225166.png', units = 'cm', width = 12, height = 10)
saveRDS(human_hca, 'tmp/integrated.human_hca.alignment.Rds')


### iter resolutions ###
dir.create('res')
dir.create('res/compare')
human_hca <- FindNeighbors(object = human_hca, dims = 1:59)
for (res in c(1:15)){
  res <- res/10
  human_hca <- FindClusters(object = human_hca, resolution = res)
  
  DimPlot(human_hca, label = T, dims = c(1,2), reduction = 'umap')
  ggsave(paste(c('res/compare/integrated.umap.1_2.res_',res,'.png'), collapse = ''), units = 'cm', width = 14, height = 10)
  
}


### fix resoultion ###
library(Seurat)
library(ggplot2)
library(plyr)
library(dplyr)

human_hca <- readRDS('tmp/integrated.human_hca.alignment.Rds')
human_hca <- FindNeighbors(object = human_hca, dims = 1:59)
human_hca <- FindClusters(object = human_hca, resolution = 0.7)

DimPlot(human_hca, reduction = "umap", pt.size = .5, label = T)
ggsave('res/integrated.umap.1_2.min_dist_0.3.res_0.7.png', units = 'cm', width = 14, height = 10)
DimPlot(human_hca, reduction = "umap", pt.size = .5, dims = c(1, 3), label = T)
ggsave('res/integrated.umap.1_3.min_dist_0.3.res_0.7.png', units = 'cm', width = 14, height = 10)
DimPlot(human_hca, reduction = "umap", pt.size = .5, dims = c(2, 3), label = T)
ggsave('res/integrated.umap.2_3.min_dist_0.3.res_0.7.png', units = 'cm', width = 14, height = 10)
saveRDS(human_hca, 'tmp/integrated.human_hca.alignment.Rds')


### re-name clusters ###
library(Seurat)
library(ggplot2)
library(plyr)
library(dplyr)


### hca ###
hca <- readRDS('tmp/integrated.human_hca.alignment.Rds')
head(hca@meta.data)
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




