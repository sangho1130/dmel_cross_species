library(Seurat)
library(plyr)

### tabula muris 10X ###
tm_tenx <- readRDS('../tabulamuris_10x/tmp/tm_10x.harmony.clustered.Rds')

### tabula muris SS2 ###
tm_ss2 <- readRDS('../tabulamuris_ss2/tmp/tm_ss2.harmony.clustered.Rds')

### mca ###
mca <- readRDS('../mca/tmp/mca.harmony.clustered.v3.noNeu.Rds')

head(tm_tenx@meta.data, n=3)
head(tm_ss2@meta.data, n=3)
head(mca@meta.data, n=3)

tm_tenx@meta.data$orig.ident <- 'TM10X'
tm_tenx@meta.data$RNA_snn_res.0.5 <- NULL
tm_tenx@meta.data$seurat_clusters <- NULL

tm_ss2@meta.data$orig.ident <- 'TMSS2'
tm_ss2@meta.data$RNA_snn_res.0.6 <- NULL
tm_ss2@meta.data$seurat_clusters <- NULL

mca@meta.data$orig.ident <- 'MCA'
mca@meta.data$RNA_snn_res.0.3 <- NULL
mca@meta.data$seurat_clusters <- NULL

head(tm_tenx@meta.data, n=3)
head(tm_ss2@meta.data, n=3)
head(mca@meta.data, n=3)

### merge ###
dir.create('tmp')

length(intersect(rownames(tm_tenx), rownames(tm_ss2))) # 15618
length(intersect(rownames(mca), intersect(rownames(tm_tenx), rownames(tm_ss2)))) # 11504
common_genes <- intersect(rownames(mca), intersect(rownames(tm_tenx), rownames(tm_ss2)))

mca_flt <- subset(mca, features = common_genes); mca_flt
mca_flt <- NormalizeData(mca_flt, normalization.method = "LogNormalize", scale.factor = 10000)

tm_tenx_flt <- subset(tm_tenx, features = common_genes); tm_tenx_flt
tm_tenx_flt <- NormalizeData(tm_tenx_flt, normalization.method = "LogNormalize", scale.factor = 10000)


tm_ss2_flt <- subset(tm_ss2, features = common_genes); tm_ss2_flt

tm_ss2_flt_mat <- as.matrix(GetAssayData(tm_ss2_flt, slot = 'data'))
summary( colSums(exp(tm_ss2_flt_mat)-1) )
factors <- 10**6/colSums(exp(tm_ss2_flt_mat)-1)
tm_ss2_flt_mat <- exp(tm_ss2_flt_mat) -1

for (n in c(1:ncol(tm_ss2_flt_mat)) ) {
  tm_ss2_flt_mat[, n] <- tm_ss2_flt_mat[, n]*factors[n]
}
summary( colSums(tm_ss2_flt_mat) ) # 10^6
tm_ss2_flt_mat <- tm_ss2_flt_mat/100
summary( colSums(tm_ss2_flt_mat) ) # 10^4
tm_ss2_flt_mat <- log(tm_ss2_flt_mat+1)
#summary( colSums(exp(tm_ss2_flt_mat)-1) )

head(tm_ss2_flt@meta.data)
tm_ss2_flt_new <- CreateSeuratObject(counts = tm_ss2_flt_mat, project = "TMSS2")
tm_ss2_flt_new <- AddMetaData(object = tm_ss2_flt_new, metadata = tm_ss2_flt@meta.data[, -c(2,3)])
head(tm_ss2_flt_new@meta.data)

mouse_tm <- merge(x = tm_tenx_flt, y = tm_ss2_flt_new)
mouse <- merge(x = mca_flt, y = mouse_tm)
head(mouse@meta.data)
remove(mouse_tm, mca, mca_flt, tm_tenx, tm_tenx_flt, tm_ss2, tm_ss2_flt, tm_ss2_flt_new, common_genes)


mouse@meta.data$orig.ident <- factor(mouse@meta.data$orig.ident, levels = c('MCA', 'TM10X', 'TMSS2'))
mouse@meta.data$mouse.id <- factor(mouse@meta.data$mouse.id, levels = unique(mouse@meta.data$mouse.id))
mouse@meta.data$dataset <- factor(mouse@meta.data$dataset, levels = unique(mouse@meta.data$dataset))
mouse@meta.data$mouse.sex <- factor(mouse@meta.data$mouse.sex, levels = unique(mouse@meta.data$mouse.sex))
mouse@meta.data$tissue <- factor(mouse@meta.data$tissue, levels = unique(mouse@meta.data$tissue))
mouse@meta.data$celltype_refined <- factor(mouse@meta.data$celltype_refined, 
                                           levels = c("Progenitors", "Granulocytopoietic cells", "Neutrophils", "Erythroids", 
                                                      "Monocytes", "Macrophage", "DC", "pDC", "Basophil",
                                                      "T cells", "NK cells", "B cells", "Plasma cells"))
head(mouse@meta.data)
summary(mouse@meta.data)
summary(mouse@meta.data$celltype_refined)
Idents(mouse) <- 'celltype_refined'; levels(Idents(mouse))
saveRDS(mouse, 'tmp/mouse.Rds')


### Standard workflow ###
library(ggplot2)
library(dplyr)
library(harmony)

dir.create('stats')
dir.create('tsne')
dir.create('tsne/compare')
dir.create('umap')
dir.create('umap/compare')
dir.create('degs')

#mouse <- readRDS('tmp/mouse.Rds')
head(mouse@meta.data); nrow(mouse@meta.data) # 17854 cells

mouse <- FindVariableFeatures(object = mouse, selection.method = "vst", nfeatures = 2000)
mouse <- ScaleData(object = mouse, vars.to.regress = c('nCount_RNA'))
mouse <- RunPCA(object = mouse, npcs = 100)
mouse <- JackStraw(object = mouse, num.replicate = 100, dims = 100)
mouse <- ScoreJackStraw(object = mouse, dims = 1:100)
JackStrawPlot(mouse, dims = 1:100) # 50 PCs 
ggsave('stats/2-4.pca.JackStrawPlot.mitoCut.pdf', units = 'cm', width = 35, height = 12)
saveRDS(mouse, 'tmp/mouse_2.Rds')

### 
summary(mouse$dataset)
mouse <- RunHarmony(mouse, group.by.vars = "dataset")
saveRDS(mouse, 'tmp/mouse.harmony.Rds')


#############
library(Seurat)
library(ggplot2)
library(plyr)
library(dplyr)
library(harmony)

#mouse <- readRDS('tmp/mouse.harmony.Rds')
### UMAP ###
for (seed in c(210728250:210728299)){
  mouse <- RunUMAP(mouse, dims=1:50, reduction = "harmony", reduction.key='UMAP', n.components=3, min.dist=0.4, seed.use = seed)
  DimPlot(mouse, group.by = "celltype_refined", reduction = "umap", label = T, label.size = 2)
  ggsave(paste0(c('umap/compare/umap.1_2.celltype_refined.seed', seed, '.png'), collapse = ''), units = 'cm', width = 16, height = 10)
} # 210728176; 210728184; 210728191; 210728268

### tSNE ###
for (seed in c(210728100:210728149)){
  mouse <- RunTSNE(mouse, dims=1:50, reduction = "harmony", reduction.key='tSNE', dim.embed=3, seed.use = seed)
  DimPlot(mouse, group.by = "celltype_refined", reduction = "tsne", label = T, label.size = 2)
  ggsave(paste0(c('tsne/compare/tsne.1_2.celltype_refined.seed', seed, '.png'), collapse = ''), units = 'cm', width = 16, height = 10)
} # 210728138


mouse <- RunUMAP(mouse, dims=1:50, reduction = "harmony", reduction.key='UMAP', n.components=3, min.dist=0.4, seed.use = 210728268)
DimPlot(object=mouse, group.by = "dataset", reduction = 'umap')
ggsave('umap/umap.1_2.dataset.seed210728268.pdf', units = 'cm', width = 13.5, height = 10)
DimPlot(object=mouse, group.by = "tissue", reduction = 'umap')
ggsave('umap/umap.1_2.tissue.seed210728268.pdf', units = 'cm', width = 14.5, height = 10)
DimPlot(object=mouse, group.by = "celltype_refined", reduction = 'umap', label = T, label.size = 2)
ggsave('umap/umap.1_2.celltype_refined.seed210728268.pdf', units = 'cm', width = 16, height = 10)
DimPlot(object=mouse, group.by = "celltype_refined", reduction = 'umap', dims = c(1,3), label = T, label.size = 2)
ggsave('umap/umap.1_3.celltype_refined.seed210728268.pdf', units = 'cm', width = 16, height = 10)
DimPlot(object=mouse, group.by = "celltype_refined", reduction = 'umap', dims = c(2,3), label = T, label.size = 2)
ggsave('umap/umap.2_3.celltype_refined.seed210728268.pdf', units = 'cm', width = 16, height = 10)

mouse <- RunTSNE(mouse, dims=1:50, reduction = "harmony", reduction.key='tSNE', dim.embed=3, seed.use = 210728138)
DimPlot(object=mouse, group.by = "dataset", reduction = 'tsne')
ggsave('tsne/tsne.1_2.dataset.seed210728138.pdf', units = 'cm', width = 13.5, height = 10)
DimPlot(object=mouse, group.by = "tissue", reduction = 'tsne')
ggsave('tsne/tsne.1_2.tissue.seed210728138.pdf', units = 'cm', width = 14.5, height = 10)
DimPlot(object=mouse, group.by = "celltype_refined", reduction = 'tsne')
ggsave('tsne/tsne.1_2.celltype_refined.seed210728138.pdf', units = 'cm', width = 16, height = 10)
saveRDS(mouse, 'tmp/mouse.harmony.Rds')


### markers ###
Idents(mouse) <- 'celltype_refined'
markers <- FindAllMarkers(object = mouse, min.pct = 0.25, logfc.threshold = 0.25, only.pos = T, test.use = 'MAST')
markers %>% group_by(cluster) %>% top_n(2, avg_logFC)
write.table(markers, 'degs/markers.mouse.celltype_refined.txt', sep = '\t', quote = F, col.names = T, row.names = F) 
top10 <- markers %>% group_by(cluster) %>% top_n(10, avg_logFC)
DoHeatmap(object = mouse, features = top10$gene, angle = 90, size = 5, raster = T, draw.lines = F)
ggsave('degs/markers.mouse.celltype_refined.pdf', units = 'cm', width = 40, height = 30)



###
library(Seurat)
library(plyr)
library(dplyr)
library(ggplot2)

mouse <- readRDS('tmp/mouse.harmony.Rds')
head(mouse@meta.data, n=3)

mouse@meta.data$atlas <- mapvalues(mouse@meta.data$orig.ident, from = levels(mouse@meta.data$orig.ident), to = c('MCA', 'TM', 'TM'))
DimPlot(object=mouse, group.by = "orig.ident", reduction = 'umap') +
  facet_wrap(~mouse@meta.data$celltype_refined, ncol = 4)
ggsave('umap/umap.1_2.dataset.by_celltype_refined.seed210728268.png', units = 'cm', width = 20, height = 20)

DimPlot(object=mouse, group.by = "orig.ident", reduction = 'umap', dims = c(1,3)) +
  facet_wrap(~mouse@meta.data$celltype_refined, ncol = 4)
ggsave('umap/umap.1_3.dataset.by_celltype_refined.seed210728268.png', units = 'cm', width = 20, height = 20)


DimPlot(subset(mouse, celltype_refined == 'Progenitors'), group.by = "orig.ident", reduction = 'umap') +
  facet_wrap(~subset(mouse, celltype_refined == 'Progenitors')@meta.data$orig.ident, ncol = 4)

DimPlot(subset(mouse, celltype_refined == 'Progenitors'), group.by = "orig.ident", reduction = 'umap', dims = c(1,3))
DimPlot(subset(mouse, celltype_refined == 'Progenitors'), group.by = "orig.ident", reduction = 'umap', dims = c(1,3)) +
  facet_wrap(~subset(mouse, celltype_refined == 'Progenitors')@meta.data$orig.ident, ncol = 4)


### Filtration ###
umapcoord <- data.frame(Embeddings(mouse, reduction = 'umap'), check.rows = F, check.names = F)
umapcoord$celltype <- mouse@meta.data$celltype_refined
head(umapcoord)

dir.create('doublets')
removecells <- c()

# Prog.
prog <- subset(umapcoord, celltype == 'Progenitors')
prog$cols <- 'singlet'
prog[rownames(subset(prog, UMAP_1 > -2 | UMAP_3 > 1.5 | UMAP_3 < -4.2)), 'cols'] <- 'doublet'

ggplot(prog, aes(UMAP_1, UMAP_3, col = cols)) + 
  geom_point() +
  scale_color_manual(values = c('red2', 'black')) +
  geom_hline(yintercept = -4.2) +
  geom_hline(yintercept = 1.5) +
  geom_vline(xintercept = -2) +
  theme_bw()
ggsave('doublets/umap.prog.png', units = 'cm', width = 11, height = 8)
removecells <- append(removecells, c( rownames(subset(prog, cols == 'doublet')) ))


# Gran. prog.
granprog <- subset(umapcoord, celltype == 'Granulocytopoietic cells')
granprog$cols <- 'singlet'
granprog[rownames(subset(granprog, UMAP_3 > -4.2 | UMAP_1 > -2.5)), 'cols'] <- 'doublet'

ggplot(granprog, aes(UMAP_1, UMAP_3, col = cols)) + 
  geom_point() +
  scale_color_manual(values = c('red2', 'black')) +
  geom_hline(yintercept = -4.2) +
  geom_vline(xintercept = -2.5) +
  theme_bw()
ggsave('doublets/umap.granprog.png', units = 'cm', width = 11, height = 8)
removecells <- append(removecells, c( rownames(subset(granprog, cols == 'doublet')) ))


# Neu
neu <- subset(umapcoord, celltype == 'Neutrophils')
neu$cols <- 'singlet'
neu[rownames(subset(neu, UMAP_3 > -5)), 'cols'] <- 'doublet'

ggplot(neu, aes(UMAP_1, UMAP_3, col = cols)) + 
  geom_point() +
  scale_color_manual(values = c('red2', 'black')) +
  geom_hline(yintercept = -5) +
  theme_bw()
ggsave('doublets/umap.neu.png', units = 'cm', width = 11, height = 8)
removecells <- append(removecells, c( rownames(subset(neu, cols == 'doublet')) ))


# Ery
ery <- subset(umapcoord, celltype == 'Erythroids')
ery$cols <- 'singlet'
ery[rownames(subset(ery, UMAP_3 < 1.2 | UMAP_1 > 3)), 'cols'] <- 'doublet'

ggplot(ery, aes(UMAP_1, UMAP_3, col = cols)) + 
  geom_point() +
  scale_color_manual(values = c('red2', 'black')) +
  geom_hline(yintercept = 1.2) +
  geom_vline(xintercept = 3) +
  theme_bw()
ggsave('doublets/umap.ery.png', units = 'cm', width = 11, height = 8)
removecells <- append(removecells, c( rownames(subset(ery, cols == 'doublet')) ))


# mono
mono <- subset(umapcoord, celltype == 'Monocytes')
mono$cols <- 'singlet'
mono[rownames(subset(mono, UMAP_3 < 1 | UMAP_1 > -3)), 'cols'] <- 'doublet'

ggplot(mono, aes(UMAP_1, UMAP_3, col = cols)) + 
  geom_point() +
  scale_color_manual(values = c('red2', 'black')) +
  geom_hline(yintercept = 1) +
  geom_vline(xintercept = -3) +
  theme_bw()
ggsave('doublets/umap.mono.png', units = 'cm', width = 11, height = 8)
removecells <- append(removecells, c( rownames(subset(mono, cols == 'doublet')) ))


# mac
mac <- subset(umapcoord, celltype == 'Macrophage')
mac$cols <- 'singlet'
mac[rownames(subset(mac, UMAP_3 < -5)), 'cols'] <- 'doublet'

ggplot(mac, aes(UMAP_1, UMAP_3, col = cols)) + 
  geom_point() +
  scale_color_manual(values = c('red2', 'black')) +
  geom_hline(yintercept = -5) +
  theme_bw()
ggsave('doublets/umap.mac.png', units = 'cm', width = 11, height = 8)
removecells <- append(removecells, c( rownames(subset(mac, cols == 'doublet')) ))


# dc
dc <- subset(umapcoord, celltype == 'DC')
dc$cols <- 'singlet'
dc[rownames(subset(dc, UMAP_1 > -3 | UMAP_3 < 2)), 'cols'] <- 'doublet'

ggplot(dc, aes(UMAP_1, UMAP_3, col = cols)) + 
  geom_point() +
  scale_color_manual(values = c('red2', 'black')) +
  geom_hline(yintercept = 2) +
  geom_vline(xintercept = -3) +
  theme_bw()
ggsave('doublets/umap.dc.png', units = 'cm', width = 11, height = 8)
removecells <- append(removecells, c( rownames(subset(dc, cols == 'doublet')) ))


# pdc
pdc <- subset(umapcoord, celltype == 'pDC')
pdc$cols <- 'singlet'
pdc[rownames(subset(pdc, UMAP_3 > 0 | UMAP_1 > -5)), 'cols'] <- 'doublet'

ggplot(pdc, aes(UMAP_1, UMAP_3, col = cols)) + 
  geom_point() +
  scale_color_manual(values = c('red2', 'black')) +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = -5) +
  theme_bw()
ggsave('doublets/umap.pdc.png', units = 'cm', width = 11, height = 8)
removecells <- append(removecells, c( rownames(subset(pdc, cols == 'doublet')) ))


# basophil
basophil <- subset(umapcoord, celltype == 'Basophil')
basophil$cols <- 'singlet'
basophil[rownames(subset(basophil, UMAP_3 < 5)), 'cols'] <- 'doublet'

ggplot(basophil, aes(UMAP_1, UMAP_3, col = cols)) + 
  geom_point() +
  scale_color_manual(values = c('red2', 'black')) +
  geom_hline(yintercept = 5) +
  theme_bw()
ggsave('doublets/umap.basophil.png', units = 'cm', width = 11, height = 8)
removecells <- append(removecells, c( rownames(subset(basophil, cols == 'doublet')) ))


# tcell
tcell <- subset(umapcoord, celltype %in% c('T cells'))
tcell$cols <- 'singlet'
tcell[rownames(subset(tcell, UMAP_1 < 0.5 | UMAP_3 < 0 | UMAP_1 > 7)), 'cols'] <- 'doublet'
tcell[rownames(subset(tcell, UMAP_1 < 1 & UMAP_3 < .5 )), 'cols'] <- 'doublet'

ggplot(tcell, aes(UMAP_1, UMAP_3, col = cols)) + 
  geom_point() +
  scale_color_manual(values = c('red2', 'black')) +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0.5) +
  geom_vline(xintercept = 7) +
  theme_bw()
ggsave('doublets/umap.tcell.png', units = 'cm', width = 11, height = 8)
removecells <- append(removecells, c( rownames(subset(tcell, cols == 'doublet')) ))

# nk
nk <- subset(umapcoord, celltype %in% c('NK cells'))
nk$cols <- 'singlet'
nk[rownames(subset(nk, UMAP_1 < 0 | UMAP_3 < 0)), 'cols'] <- 'doublet'

ggplot(nk, aes(UMAP_1, UMAP_3, col = cols)) + 
  geom_point() +
  scale_color_manual(values = c('red2', 'black')) +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  theme_bw()
ggsave('doublets/umap.nk.png', units = 'cm', width = 11, height = 8)
removecells <- append(removecells, c( rownames(subset(nk, cols == 'doublet')) ))


# b cell
bcell <- subset(umapcoord, celltype %in% c('B cells'))
bcell$cols <- 'singlet'
bcell[rownames(subset(bcell, UMAP_1 < 4)), 'cols'] <- 'doublet'
bcell[rownames(subset(bcell, UMAP_1 < 6 & UMAP_3 < 3 )), 'cols'] <- 'doublet'

ggplot(bcell, aes(UMAP_1, UMAP_3, col = cols)) + 
  geom_point() +
  scale_color_manual(values = c('red2', 'black')) +
  geom_vline(xintercept = 4) +
  theme_bw()
ggsave('doublets/umap.bcell.png', units = 'cm', width = 11, height = 8)
removecells <- append(removecells, c( rownames(subset(bcell, cols == 'doublet')) ))


# pc
plasma <- subset(umapcoord, celltype %in% c('Plasma cells'))
plasma$cols <- 'singlet'
plasma[rownames(subset(plasma, UMAP_1 < 2 | UMAP_3 > -3)), 'cols'] <- 'doublet'

ggplot(plasma, aes(UMAP_1, UMAP_3, col = cols)) + 
  geom_point() +
  scale_color_manual(values = c('red2', 'black')) +
  geom_vline(xintercept = 2) +
  geom_hline(yintercept = -3) +
  theme_bw()
ggsave('doublets/umap.plasma.png', units = 'cm', width = 11, height = 8)
removecells <- append(removecells, c( rownames(subset(plasma, cols == 'doublet')) ))


##
length(removecells) # 2686

umapcoord$remove <- 'no'
umapcoord[removecells, 'remove'] <- 'yes'

ggplot(umapcoord, aes(UMAP_1, UMAP_3, col = remove)) + 
  geom_point(size = .5) +
  scale_color_manual(values = c('black', 'red2')) +
  theme_bw()
ggsave('doublets/umap.allcells.png', units = 'cm', width = 10, height = 8)

ggplot(umapcoord, aes(UMAP_1, UMAP_3, col = remove)) + 
  geom_point(size = .5) +
  facet_wrap(~remove) +
  scale_color_manual(values = c('black', 'red2')) +
  theme_bw() +
  theme(legend.position = 'none')
ggsave('doublets/umap.allcells.byremove.png', units = 'cm', width = 15, height = 8)
#saveRDS(removecells, 'tmp/doublets.Rds')


mouse_flt <- subset(mouse, cells = setdiff(rownames(mouse@meta.data), removecells))
mouse_flt@meta.data <- droplevels(mouse_flt@meta.data)
head(mouse_flt@meta.data, n=3); summary(mouse_flt@meta.data)

Idents(mouse_flt) <- 'celltype_refined'
summary(Idents(mouse_flt))
#saveRDS(mouse_flt, 'tmp/mouse.harmony.flt.Rds')


DimPlot(object=mouse_flt, group.by = "dataset", reduction = 'umap')
ggsave('umap/umap_flt.1_2.dataset.seed210728268.pdf', units = 'cm', width = 13.5, height = 10)
DimPlot(object=mouse_flt, group.by = "tissue", reduction = 'umap')
ggsave('umap/umap_flt.1_2.tissue.seed210728268.pdf', units = 'cm', width = 14.5, height = 10)
DimPlot(object=mouse_flt, group.by = "celltype_refined", reduction = 'umap')
ggsave('umap/umap_flt.1_2.celltype_refined.seed210728268.pdf', units = 'cm', width = 16, height = 10)
DimPlot(object=mouse_flt, group.by = "celltype_refined", reduction = 'umap', dim = c(1,3))
ggsave('umap/umap_flt.1_3.celltype_refined.seed210728268.pdf', units = 'cm', width = 16, height = 10)
DimPlot(object=mouse_flt, group.by = "celltype_refined", reduction = 'umap', dim = c(2,3))
ggsave('umap/umap_flt.2_3.celltype_refined.seed210728268.pdf', units = 'cm', width = 16, height = 10)
DimPlot(object=mouse_flt, group.by = "orig.ident", reduction = 'umap') +
  facet_wrap(~mouse_flt@meta.data$celltype_refined, ncol = 4)
ggsave('umap/umap_flt.1_2.dataset.by_celltype_refined.seed210728268.png', units = 'cm', width = 20, height = 20)
DimPlot(object=mouse_flt, group.by = "orig.ident", reduction = 'umap', dim = c(1,3)) +
  facet_wrap(~mouse_flt@meta.data$celltype_refined, ncol = 4)
ggsave('umap/umap_flt.1_3.dataset.by_celltype_refined.seed210728268.png', units = 'cm', width = 20, height = 20)
DimPlot(object=mouse_flt, group.by = "orig.ident", reduction = 'umap', dim = c(2,3)) +
  facet_wrap(~mouse_flt@meta.data$celltype_refined, ncol = 4)
ggsave('umap/umap_flt.2_3.dataset.by_celltype_refined.seed210728268.png', units = 'cm', width = 20, height = 20)

DimPlot(object=mouse_flt, group.by = "dataset", reduction = 'tsne')
ggsave('tsne/tsne_flt.1_2.dataset.seed210728138.pdf', units = 'cm', width = 13.5, height = 10)
DimPlot(object=mouse_flt, group.by = "tissue", reduction = 'tsne')
ggsave('tsne/tsne_flt.1_2.tissue.seed210728138.pdf', units = 'cm', width = 14.5, height = 10)
DimPlot(object=mouse_flt, group.by = "celltype_refined", reduction = 'tsne')
ggsave('tsne/tsne_flt.1_2.celltype_refined.seed210728138.pdf', units = 'cm', width = 16, height = 10)


### markers ###
markers <- FindAllMarkers(object = mouse_flt, min.pct = 0.25, logfc.threshold = 0.25, only.pos = T, test.use = 'MAST')
markers %>% group_by(cluster) %>% top_n(2, avg_logFC)
write.table(markers, 'degs/markers.mouse.celltype_refined.DBflt.txt', sep = '\t', quote = F, col.names = T, row.names = F) 
top10 <- markers %>% group_by(cluster) %>% top_n(10, avg_logFC)
DoHeatmap(object = mouse_flt, features = top10$gene, angle = 90, size = 5, raster = T, draw.lines = F)
ggsave('degs/markers.mouse.celltype_refined.DBflt.pdf', units = 'cm', width = 40, height = 30)


### Figures ###
library(Seurat)
library(plyr)
library(dplyr)
library(ggplot2)
library(RColorBrewer)

mouse_flt <- readRDS('tmp/mouse.harmony.flt.Rds')
head(mouse_flt@meta.data, n=3); nrow(mouse_flt@meta.data) # 15168 cells
summary(mouse_flt@meta.data$orig.ident)


use_cols <- colorRampPalette(brewer.pal(8, "Accent"))(18)
use_cols <- c('Progenitors' = use_cols[c(1)],
              'Granulocytopoietic cells' = use_cols[c(2)],
              'Erythroids' = use_cols[c(3)], 
              'Neutrophils' = use_cols[c(4)], 
              'Monocytes' = use_cols[c(5)], 
              'Macrophage' = use_cols[c(6)], 
              'DC' = use_cols[c(7)], 
              'pDC' = use_cols[c(8)], 
              'Basophil' = use_cols[c(9)], 
              'T cells' = use_cols[c(10)], 
              'NK cells' = use_cols[c(12)], 
              'B cells' = use_cols[c(16)], 
              'Plasma cells' = use_cols[c(17)])

DimPlot(object=mouse_flt, group.by = "celltype_refined", reduction = 'umap', pt.size = .5) +
  scale_colour_manual(name="Cell type", values = use_cols) + 
  theme_void() +
  theme(legend.position = 'none')
#ggsave('umap/umap_flt.1_2.celltype_refined.seed210728268.void.png', units = 'cm', width = 6, height = 6)
DimPlot(object=mouse_flt, group.by = "celltype_refined", reduction = 'umap') +
  scale_colour_manual(name="Cell type", values = use_cols)
#ggsave('umap/umap_flt.1_2.celltype_refined.seed210728268.pdf', units = 'cm', width = 14, height = 8)




