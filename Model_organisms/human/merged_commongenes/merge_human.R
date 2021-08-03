library(Seurat)
library(plyr)

### HCA ###
hca <- readRDS('../hca/tmp/integrated.hca.alignment.filtered.Rds')

### HCL ###
hcl <- readRDS('../hcl/tmp/human_hcl.harmony.Rds')

head(hca@meta.data, n=3) # 243398
head(hcl@meta.data, n=3) # 21568


hca@meta.data$integrated_snn_res.0.7 <- NULL
hca@meta.data$seurat_clusters <- NULL
hca@meta.data$sample <- 'AdultBoneMarrow'
hca@meta.data$percent.mt <- NULL

hcl@meta.data$RNA_snn_res.0.3 <- NULL
hcl@meta.data$cluster <- NULL
colnames(hcl@meta.data)[7] <- 'Library'
hcl@meta.data$sex <- 'NS'
hcl@meta.data$celltype <- NULL
hcl@meta.data$celltype_broad <- NULL
hcl@meta.data$stage <- NULL
hcl@meta.data$batch <- NULL


hca@meta.data <- hca@meta.data[, c("orig.ident", "nCount_RNA", "nFeature_RNA", "Library", "sample", "sex", "celltype_refined")]
head(hca@meta.data, n=3) # 243398
hcl@meta.data <- hcl@meta.data[, c("orig.ident", "nCount_RNA", "nFeature_RNA", "Library", "sample", "sex", "celltype_refined")]
head(hcl@meta.data, n=3) # 21568


### merge ###
dir.create('tmp')

length(rownames(hca)) # 58347
length(rownames(hcl)) # 27341
length(intersect(rownames(hca), rownames(hcl))) # 20935
common_genes <- intersect(rownames(hca), rownames(hcl))


hca_flt <- subset(hca, features = common_genes); hca_flt
hca_flt <- NormalizeData(hca_flt, normalization.method = "LogNormalize", scale.factor = 10000)

hcl_flt <- subset(hcl, features = common_genes); hcl_flt
hcl_flt <- NormalizeData(hcl_flt, normalization.method = "LogNormalize", scale.factor = 10000)

human <- merge(x = hca_flt, y = hcl_flt)
head(human@meta.data, n=3); nrow(human@meta.data) # 264966

remove(hca, hca_flt, hcl, hcl_flt, common_genes)

human@meta.data$orig.ident <- factor(human@meta.data$orig.ident, levels = c('HCA', 'HCL'))
human@meta.data$Library <- factor(human@meta.data$Library, levels = unique(human@meta.data$Library))
human@meta.data$sample <- factor(human@meta.data$sample, levels = unique(human@meta.data$sample))
human@meta.data$sex <- factor(human@meta.data$sex, levels = unique(human@meta.data$sex))


human@meta.data$celltype_refined <- mapvalues(human@meta.data$celltype_refined,
                                              from = unique(human@meta.data$celltype_refined),
                                              to = c("Monocytes", "pre-B cells", "T cells", "NK cells", "DC", "Progenitors", "Erythroids", "B cells", "pro-B cells", 
                                                     "CD8 T cells", "pDC", "Plasma cells", "Granulocyte progenitor", "pre-PC", "Platelet", "Plasma cells", "Neutrophils", "Macrophage"))
human@meta.data$celltype_refined <- factor(human@meta.data$celltype_refined, 
                                           levels = c("Progenitors", "Granulocyte progenitor", "Erythroids", "Neutrophils", 
                                                      "Monocytes", "Macrophage", "DC", "pDC", "T cells", "CD8 T cells", "NK cells", 
                                                      "pre-PC", "pro-B cells", "pre-B cells", "B cells", "Plasma cells", "Plasma cell", "Platelet"))
head(human@meta.data)
summary(human@meta.data)
Idents(human) <- 'celltype_refined'; levels(Idents(human))
saveRDS(human, 'tmp/human.Rds')


### Standard workflow ###
library(Seurat)
library(plyr)
library(dplyr)
library(ggplot2)
library(harmony)

dir.create('stats')
dir.create('umap')
dir.create('umap/compare')
dir.create('degs')

#human <- readRDS('tmp/human.Rds')
#head(human@meta.data); nrow(human@meta.data) # 264966 cells

human <- FindVariableFeatures(object = human, selection.method = "vst", nfeatures = 2000)
human <- ScaleData(object = human, vars.to.regress = c('nCount_RNA'))
human <- RunPCA(object = human, npcs = 100)
human <- JackStraw(object = human, num.replicate = 100, dims = 100)
human <- ScoreJackStraw(object = human, dims = 1:100)
JackStrawPlot(human, dims = 1:100) # 68 PCs 
ggsave('stats/2-4.pca.JackStrawPlot.mitoCut.pdf', units = 'cm', width = 30, height = 12)
saveRDS(human, 'tmp/human_2.Rds')


human <- readRDS('tmp/human_2.Rds')

#human_dataset <- RunHarmony(human, group.by.vars = "orig.ident") # Quick-TRANSfer stage steps exceeded maximum (= 13248300)
#saveRDS(human_dataset, 'tmp/human.harmony_orig.ident.Rds')

human_library <- RunHarmony(human, group.by.vars = "Library") # Quick-TRANSfer stage steps exceeded maximum (= 13248300)
saveRDS(human_library, 'tmp/human.harmony_Library.Rds')

human <- human_library
remove(human_dataset, human_library) ###

### UMAP ###
for (seed in c(210721259:210721299)){
  human <- RunUMAP(human, dims=1:68, reduction = "harmony", reduction.key='UMAP', n.components=3, min.dist=0.25, seed.use = seed)
  DimPlot(human, group.by = "celltype_refined", reduction = "umap", label = T, label.size = 2)
  ggsave(paste0(c('umap/compare/umap.1_2.seed', seed, '.png'), collapse = ''), units = 'cm', width = 16, height = 10)
  DimPlot(human, group.by = "orig.ident", reduction = "umap")
  ggsave(paste0(c('umap/compare/umap.1_2.dataset.seed', seed, '.png'), collapse = ''), units = 'cm', width = 12, height = 10)
} # 210721269


human <- RunUMAP(human, dims=1:68, reduction = "harmony", reduction.key='UMAP', n.components=3, min.dist=0.25, seed.use = 210721269)
DimPlot(object=human, group.by = "orig.ident", reduction = 'umap')
ggsave('umap/umap.1_2.dataset.seed210721269.pdf', units = 'cm', width = 12, height = 10)
DimPlot(object=human, group.by = "sample", reduction = 'umap')
ggsave('umap/umap.1_2.sample.seed210721269.pdf', units = 'cm', width = 14.5, height = 10)
DimPlot(object=human, group.by = "celltype_refined", reduction = 'umap')
ggsave('umap/umap.1_2.celltype_refined.seed210721269.pdf', units = 'cm', width = 16, height = 10)

DimPlot(object=human, group.by = "orig.ident", reduction = 'umap') +
  facet_wrap(~human@meta.data$celltype_refined, ncol = 4)
ggsave('umap/umap.1_2.dataset.by_celltype_refined.seed210721269.png', units = 'cm', width = 30, height = 35)

saveRDS(human, 'tmp/human.harmony_Library.Rds')



### markers ###
human@meta.data <- droplevels(human@meta.data)
Idents(human) <- 'celltype_refined'
levels(Idents(human))

markers <- FindAllMarkers(object = human, min.pct = 0.25, logfc.threshold = 0.25, only.pos = T, test.use = 'MAST')
markers %>% group_by(cluster) %>% top_n(2, avg_logFC)
write.table(markers, 'degs/markers.human.celltype_refined.txt', sep = '\t', quote = F, col.names = T, row.names = F) 
top10 <- markers %>% group_by(cluster) %>% top_n(10, avg_logFC)
DoHeatmap(object = human, features = top10$gene, angle = 90, size = 5, raster = T, draw.lines = F)
ggsave('degs/markers.human.celltype_refined.png', units = 'cm', width =50, height = 40)


##################
### filtration ###
umapcoord <- data.frame(Embeddings(human, reduction = 'umap'), check.rows = F, check.names = F)
umapcoord$celltype <- human@meta.data$celltype_refined
head(umapcoord)

dir.create('doublets')
removecells <- c()


# Gran. prog.
granprog <- subset(umapcoord, celltype == 'Granulocyte progenitor')
granprog$cols <- 'singlet'
granprog[rownames(subset(granprog, UMAP_1 > 0 | UMAP_2 < -5)), 'cols'] <- 'doublet'

ggplot(granprog, aes(UMAP_1, UMAP_2, col = cols)) + 
  geom_point() +
  scale_color_manual(values = c('red2', 'black')) +
  geom_hline(yintercept = -5) +
  geom_vline(xintercept = 0) +
  theme_bw()
ggsave('doublets/umap.granprog.png', units = 'cm', width = 11, height = 8)
removecells <- append(removecells, c( rownames(subset(granprog, UMAP_1 > 0 | UMAP_2 < -5)) ))


# Erythroids
ery <- subset(umapcoord, celltype == 'Erythroids')
ery$cols <- 'singlet'
ery[rownames(subset(ery, UMAP_1 > 0 | UMAP_2 < -5)), 'cols'] <- 'doublet'

ggplot(ery, aes(UMAP_1, UMAP_2, col = cols)) + 
  geom_point() +
  scale_color_manual(values = c('red2', 'black')) +
  geom_hline(yintercept = -5) +
  geom_vline(xintercept = 0) +
  theme_bw()
ggsave('doublets/umap.ery.png', units = 'cm', width = 11, height = 8)
removecells <- append(removecells, c( rownames(subset(ery, UMAP_1 > 0 | UMAP_2 < -5)) ))


# Neu.
neu <- subset(umapcoord, celltype == 'Neutrophils')
neu$cols <- 'singlet'
neu[rownames(subset(neu, UMAP_1 > 0 | UMAP_2 < -5)), 'cols'] <- 'doublet'

ggplot(neu, aes(UMAP_1, UMAP_2, col = cols)) + 
  geom_point() +
  scale_color_manual(values = c('red2', 'black')) +
  geom_hline(yintercept = -5) +
  geom_vline(xintercept = 0) +
  theme_bw()
ggsave('doublets/umap.neu.png', units = 'cm', width = 11, height = 8)
removecells <- append(removecells, c( rownames(subset(neu, UMAP_1 > 0 | UMAP_2 < -5)) ))


# Mono
mono <- subset(umapcoord, celltype == 'Monocytes')
mono$cols <- 'singlet'
mono[rownames(subset(mono, UMAP_1 > 0 | UMAP_2 < 0.5)), 'cols'] <- 'doublet'

ggplot(mono, aes(UMAP_1, UMAP_2, col = cols)) + 
  geom_point() +
  scale_color_manual(values = c('red2', 'black')) +
  geom_hline(yintercept = 0.5) +
  geom_vline(xintercept = 0) +
  theme_bw()
ggsave('doublets/umap.mono.png', units = 'cm', width = 11, height = 8)
removecells <- append(removecells, c( rownames(subset(mono, UMAP_1 > 0 | UMAP_2 < 0.5)) ))


# Mac
mac <- subset(umapcoord, celltype == 'Macrophage')
mac$cols <- 'singlet'
mac[rownames(subset(mac, UMAP_1 > 0 | UMAP_2 < 0.5)), 'cols'] <- 'doublet'

ggplot(mac, aes(UMAP_1, UMAP_2, col = cols)) + 
  geom_point() +
  scale_color_manual(values = c('red2', 'black')) +
  geom_hline(yintercept = 0.5) +
  geom_vline(xintercept = 0) +
  theme_bw()
ggsave('doublets/umap.mac.png', units = 'cm', width = 11, height = 8)
removecells <- append(removecells, c( rownames(subset(mac, UMAP_1 > 0 | UMAP_2 < 0.5)) ))


# DC
dc <- subset(umapcoord, celltype == 'DC')
dc$cols <- 'singlet'
dc[rownames(subset(dc, UMAP_2 < 0.5)), 'cols'] <- 'doublet'

ggplot(dc, aes(UMAP_1, UMAP_2, col = cols)) + 
  geom_point() +
  scale_color_manual(values = c('red2', 'black')) +
  geom_hline(yintercept = 0.5) +
  theme_bw()
ggsave('doublets/umap.dc.png', units = 'cm', width = 11, height = 8)
removecells <- append(removecells, c( rownames(subset(dc, UMAP_2 < 0.5)) ))


# pdc
pdc <- subset(umapcoord, celltype == 'pDC')
pdc$cols <- 'singlet'
pdc[rownames(subset(pdc, UMAP_1 > 0 | UMAP_2 < -0.5)), 'cols'] <- 'doublet'

ggplot(pdc, aes(UMAP_1, UMAP_2, col = cols)) + 
  geom_point() +
  scale_color_manual(values = c('red2', 'black')) +
  geom_hline(yintercept = -0.5) +
  geom_vline(xintercept = 0) +
  theme_bw()
ggsave('doublets/umap.pdc.png', units = 'cm', width = 11, height = 8)
removecells <- append(removecells, c( rownames(subset(pdc, UMAP_1 > 0 | UMAP_2 < -0.5)) ))


# T/CD8/NK cells
tnk <- subset(umapcoord, celltype %in% c('T cells', 'CD8 T cells', 'NK cells'))
tnk$cols <- 'singlet'
tnk[rownames(subset(tnk, UMAP_1 < 0)), 'cols'] <- 'doublet'

ggplot(tnk, aes(UMAP_1, UMAP_2, col = cols)) + 
  geom_point() +
  scale_color_manual(values = c('red2', 'black')) +
  geom_vline(xintercept = 0) +
  theme_bw()
ggsave('doublets/umap.tnk.png', units = 'cm', width = 11, height = 8)
removecells <- append(removecells, c( rownames(subset(tnk, UMAP_1 < 0)) ))


# pre-PC
ppc <- subset(umapcoord, celltype == 'pre-PC')
ppc$cols <- 'singlet'
ppc[rownames(subset(ppc, UMAP_2 > 1 | UMAP_2 < -5)), 'cols'] <- 'doublet'

ggplot(ppc, aes(UMAP_1, UMAP_2, col = cols)) + 
  geom_point() +
  scale_color_manual(values = c('red2', 'black')) +
  geom_hline(yintercept = 1) +
  geom_hline(yintercept = -5) +
  theme_bw()
ggsave('doublets/umap.ppc.png', units = 'cm', width = 11, height = 8)
removecells <- append(removecells, c( rownames(subset(ppc, UMAP_2 > 1 | UMAP_2 < -5)) ))


# pro-B
prob <- subset(umapcoord, celltype == 'pro-B cells')
prob$cols <- 'singlet'
prob[rownames(subset(prob, UMAP_2 < -5)), 'cols'] <- 'doublet'

ggplot(prob, aes(UMAP_1, UMAP_2, col = cols)) + 
  geom_point() +
  scale_color_manual(values = c('red2', 'black')) +
  geom_hline(yintercept = -5) +
  theme_bw()
ggsave('doublets/umap.prob.png', units = 'cm', width = 11, height = 8)
removecells <- append(removecells, c( rownames(subset(prob, UMAP_2 < -5)) ))


# pre-B
preb <- subset(umapcoord, celltype == 'pre-B cells')
preb$cols <- 'singlet'
preb[rownames(subset(preb, UMAP_2 > -5)), 'cols'] <- 'doublet'

ggplot(preb, aes(UMAP_1, UMAP_2, col = cols)) + 
  geom_point() +
  scale_color_manual(values = c('red2', 'black')) +
  geom_hline(yintercept = -5) +
  theme_bw()
ggsave('doublets/umap.preb.png', units = 'cm', width = 11, height = 8)
removecells <- append(removecells, c( rownames(subset(preb, UMAP_2 > -5)) ))


# B
bcell <- subset(umapcoord, celltype == 'B cells')
bcell$cols <- 'singlet'
bcell[rownames(subset(bcell, UMAP_2 > -7)), 'cols'] <- 'doublet'

ggplot(bcell, aes(UMAP_1, UMAP_2, col = cols)) + 
  geom_point() +
  scale_color_manual(values = c('red2', 'black')) +
  geom_hline(yintercept = -7) +
  theme_bw()
ggsave('doublets/umap.bcell.png', units = 'cm', width = 11, height = 8)
removecells <- append(removecells, c( rownames(subset(bcell, UMAP_2 > -7)) ))


# Plasma cell
pc <- subset(umapcoord, celltype == 'Plasma cells')
pc$cols <- 'singlet'
pc[rownames(subset(pc, UMAP_2 > -3 | UMAP_2 < -7)), 'cols'] <- 'doublet'

ggplot(pc, aes(UMAP_1, UMAP_2, col = cols)) + 
  geom_point() +
  scale_color_manual(values = c('red2', 'black')) +
  geom_hline(yintercept = -7) +
  geom_hline(yintercept = -3) +
  theme_bw()
ggsave('doublets/umap.pc.png', units = 'cm', width = 11, height = 8)
removecells <- append(removecells, c( rownames(subset(pc, UMAP_2 > -3 | UMAP_2 < -7)) ))


# Platelet
plat <- subset(umapcoord, celltype == 'Platelet')
plat$cols <- 'singlet'
plat[rownames(subset(plat, UMAP_2 < 0.55 | UMAP_1 < 0)), 'cols'] <- 'doublet'

ggplot(plat, aes(UMAP_1, UMAP_2, col = cols)) + 
  geom_point() +
  scale_color_manual(values = c('red2', 'black')) +
  geom_hline(yintercept = 0.5) +
  geom_vline(xintercept = 0) +
  theme_bw()
ggsave('doublets/umap.plat.png', units = 'cm', width = 11, height = 8)
removecells <- append(removecells, c( rownames(subset(plat, UMAP_2 < 0.55 | UMAP_1 < 0)) ))


##
length(removecells) # 2146

umapcoord$remove <- 'no'
umapcoord[removecells, 'remove'] <- 'yes'

ggplot(umapcoord, aes(UMAP_1, UMAP_2, col = remove)) + 
  geom_point(size = .5) +
  scale_color_manual(values = c('black', 'red2')) +
  theme_bw()
ggsave('doublets/umap.allcells.png', units = 'cm', width = 10, height = 8)

#saveRDS(removecells, 'tmp/doublets.Rds')


human_flt <- subset(human, cells = setdiff(rownames(human@meta.data), removecells))
human_flt@meta.data <- droplevels(human_flt@meta.data)
human_flt@meta.data$nCount_integrated <- NULL
human_flt@meta.data$nFeature_integrated <- NULL
head(human_flt@meta.data, n=3)
saveRDS(human_flt, 'tmp/human.harmony_Library.flt.Rds')


library(RColorBrewer)
mycol <- colorRampPalette(brewer.pal(9, "Spectral"))(17)

DimPlot(human_flt, pt.size = 0.5, reduction = 'umap') + 
  scale_color_manual(values = mycol ) + 
  theme_void() + 
  theme(legend.position = 'none')
ggsave('umap/umap_flt.1_2.celltype_refined.seed210721269.void.png', units = 'cm', width = 8, height = 8)
DimPlot(human_flt, pt.size = .5, reduction = 'umap', label = T, label.size = 2) + 
  scale_color_manual(values = mycol )
ggsave('umap/umap_flt.1_2.celltype_refined.seed210721269.png', units = 'cm', width = 16, height = 10)


DimPlot(object=human_flt, group.by = "orig.ident", reduction = 'umap')
ggsave('umap/umap_flt.1_2.dataset.seed210721269.pdf', units = 'cm', width = 12, height = 10)
DimPlot(object=human_flt, group.by = "sample", reduction = 'umap')
ggsave('umap/umap_flt.1_2.sample.seed210721269.pdf', units = 'cm', width = 14.5, height = 10)

DimPlot(object=human_flt, group.by = "orig.ident", reduction = 'umap') +
  facet_wrap(~human_flt@meta.data$celltype_refined, ncol = 4)
ggsave('umap/umap_flt.1_2.dataset.by_celltype_refined.seed210721269.png', units = 'cm', width = 30, height = 35)



### markers ###
Idents(human_flt) <- 'celltype_refined'
levels(Idents(human_flt))

markers <- FindAllMarkers(object = human_flt, min.pct = 0.25, logfc.threshold = 0.25, only.pos = T, test.use = 'MAST')
markers %>% group_by(cluster) %>% top_n(2, avg_logFC)
write.table(markers, 'degs/markers.human.celltype_refined.DBflt.txt', sep = '\t', quote = F, col.names = T, row.names = F) 
top10 <- markers %>% group_by(cluster) %>% top_n(10, avg_logFC)
DoHeatmap(object = human_flt, features = top10$gene, angle = 90, size = 5, raster = T, draw.lines = F)
ggsave('degs/markers.human.celltype_refined.DBflt.png', units = 'cm', width =50, height = 40)


library(Seurat)
library(plyr)
library(dplyr)
library(ggplot2)

human_flt <- readRDS('tmp/human.harmony_Library.flt.Rds')
human_flt <- subset(human_flt, celltype_refined != "Platelet")
human_flt@meta.data <- droplevels(human_flt@meta.data)
summary(Idents(human_flt))

markers <- FindAllMarkers(object = human_flt, min.pct = 0.25, logfc.threshold = 0.25, only.pos = T, test.use = 'MAST')
markers %>% group_by(cluster) %>% top_n(2, avg_logFC)
write.table(markers, 'degs/markers.human.celltype_refined.DBflt.noPlt.txt', sep = '\t', quote = F, col.names = T, row.names = F) 
top10 <- markers %>% group_by(cluster) %>% top_n(10, avg_logFC)
DoHeatmap(object = human_flt, features = top10$gene, angle = 90, size = 5, raster = T, draw.lines = F)
ggsave('degs/markers.human.celltype_refined.DBflt.noPlt.png', units = 'cm', width =50, height = 40)


summary(human_flt@meta.data$celltype_refined)
human_flt@meta.data$celltype_refined_broad <- mapvalues(human_flt@meta.data$celltype_refined, from = levels(human_flt@meta.data$celltype_refined),
                                                        to = c( "Progenitors", "Granulocyte progenitor", "Erythroids", "Neutrophils", "MonoMac", "MonoMac", "DC", "pDC", 
                                                                "T cells", "T cells", "NK cells", "pre-PC", "pro-B cells", "pre-B cells", "B cells", "Plasma cells"))
Idents(human_flt) <- 'celltype_refined_broad'
markers <- FindAllMarkers(object = human_flt, min.pct = 0.25, logfc.threshold = 0.25, only.pos = T, test.use = 'MAST')
markers %>% group_by(cluster) %>% top_n(2, avg_logFC)
write.table(markers, 'degs/markers.human.celltype_refined_broad.DBflt.noPlt.txt', sep = '\t', quote = F, col.names = T, row.names = F) 
top10 <- markers %>% group_by(cluster) %>% top_n(10, avg_logFC)
DoHeatmap(object = human_flt, features = top10$gene, angle = 90, size = 5, raster = T, draw.lines = F)
ggsave('degs/markers.human.celltype_refined_broad.DBflt.noPlt.pdf', units = 'cm', width =50, height = 40)


### Figures ###
library(Seurat)
library(plyr)
library(dplyr)
library(ggplot2)
library(RColorBrewer)

human_flt <- readRDS('tmp/human.harmony_Library.flt.Rds')
human_flt <- subset(human_flt, celltype_refined != "Platelet")
human_flt@meta.data <- droplevels(human_flt@meta.data)
summary(Idents(human_flt))
head(human_flt@meta.data, n=3); nrow(human_flt@meta.data) # 262630 cells

use_cols <- colorRampPalette(brewer.pal(8, "Accent"))(18)
use_cols <- c('Progenitors' = use_cols[c(1)],
              'Granulocyte progenitor' = use_cols[c(2)],
              'Erythroids' = use_cols[c(3)],
              'Neutrophils' = use_cols[c(4)],
              'Monocytes' = use_cols[c(5)], 
              'Macrophage' = use_cols[c(6)], 
              'DC' = use_cols[c(7)], 
              'pDC' = use_cols[c(8)], 
              'T cells' = use_cols[c(10)], 
              'CD8 T cells' = use_cols[c(11)], 
              'NK cells' = use_cols[c(12)], 
              'pre-PC' = use_cols[c(13)],
              'pro-B cells' = use_cols[c(14)],
              'pre-B cells' = use_cols[c(15)],
              'B cells' = use_cols[c(16)], 
              'Plasma cells' = use_cols[c(17)])

DimPlot(object=human_flt, group.by = "celltype_refined", reduction = 'umap', pt.size = .5) +
  scale_colour_manual(name="Cell type", values = use_cols) + 
  theme_void() +
  theme(legend.position = 'none')
ggsave('umap/umap_flt.1_2.celltype_refined.seed210721269.void.png', units = 'cm', width = 6, height = 6)
DimPlot(object=human_flt, group.by = "celltype_refined", reduction = 'umap') +
  scale_colour_manual(name="Cell type", values = use_cols)
ggsave('umap/umap_flt.1_2.celltype_refined.seed210721269.png', units = 'cm', width = 14, height = 8)


DimPlot(object=human_flt, group.by = "celltype_refined", reduction = 'umap', pt.size = .5, dims = c(1,3)) +
  scale_colour_manual(name="Cell type", values = use_cols) + 
  theme_void() +
  theme(legend.position = 'none')
ggsave('umap/umap_flt.1_3.celltype_refined.seed210721269.void.png', units = 'cm', width = 6, height = 6)

DimPlot(object=human_flt, group.by = "celltype_refined", reduction = 'umap', pt.size = .5, dims = c(2,3)) +
  scale_colour_manual(name="Cell type", values = use_cols) + 
  theme_void() +
  theme(legend.position = 'none')
ggsave('umap/umap_flt.2_3.celltype_refined.seed210721269.void.png', units = 'cm', width = 6, height = 6)
