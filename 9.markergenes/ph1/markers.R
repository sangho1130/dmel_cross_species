
library(Seurat)
library(plyr)
library(ggplot2)

fly <- readRDS('../../../Drop-seq_merge_harmony_6.22/tmp/dropseq.combined.harmony_flt.Rds')
Idents(fly) <- 'celltype'
bloodcells <- subset(fly, celltype != 'PSC'); bloodcells@meta.data <- droplevels(bloodcells@meta.data)
remove(fly)
bloodcells@meta.data$celltype_v2 <- mapvalues(bloodcells@meta.data$subclustering,
                                              from = levels(bloodcells@meta.data$subclustering),
                                              to = c("PH 1", "PH", "PH", "PH", "PH", "PH", "PM", "PM 120", "PM 120", "PM 120", 
                                                     "LM", "LM", "CC", "CC", "GST-rich", "Adipohemocyte"))
Idents(bloodcells) <- 'celltype_v2'
levels(Idents(bloodcells))


fish <- readRDS('../../zebrafish/Tang_Q_et_al/InDrops/tmp/zebrafish.harmony.bloodcells.Rds')
levels(Idents(fish))

mca <- readRDS('../../mouse/mca/tmp/mca.harmony.clustered.v3.noNeu.Rds')
levels(Idents(mca))

tm_merged <- readRDS('../../mouse/tabulamuris_merged/tmp/tm_merged.harmony.Rds')
Idents(tm_merged) <- 'celltype_refined'

tm_10x <- readRDS('../../mouse/tabulamuris_10x/tmp/tm_10x.harmony.clustered.Rds')
levels(Idents(tm_10x))
tm_ss2 <- readRDS('../../mouse/tabulamuris_ss2/tmp/tm_ss2.harmony.clustered.Rds')
levels(Idents(tm_ss2))

hcl <- readRDS('../../human/hcl/tmp/human_hcl.harmony.Rds')
levels(Idents(hcl))


###
flygene <- c('Dsp1', 'smt3', 'Set', 'PCNA', 'Fs(2)Ket', 'hth')
fishgene <- c('hmgb2b', 'sumo3a', 'seta', 'pcna', 'kpnb1', 'meis1b')
mousegene <- c('Hmgb3', 'Sumo3', 'Set', 'Pcna', 'Kpnb1', 'Meis1')
humangene <- c('HMGB3', 'SUMO3', 'SET', 'PCNA', 'KPNB1', 'MEIS1')


DotPlot(bloodcells, features = flygene) + labs(x = '', y = '') + coord_flip() + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5))
#ggsave('drosophila.pdf', units = 'cm', width = 11, height = 10)

DotPlot(fish, features = fishgene) + labs(x = '', y = '')  + coord_flip() + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5))
#ggsave('zebrafish.pdf', units = 'cm', width = 11, height = 10)

DotPlot(mca, features = mousegene) + labs(x = '', y = '')  + coord_flip() + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5))
#ggsave('mouse_mca.pdf', units = 'cm', width = 14, height = 10)

DotPlot(tm_merged, features = mousegene) + labs(x = '', y = '')  + coord_flip() + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5))
#ggsave('mouse_tm.pdf', units = 'cm', width = 14, height = 10)
DotPlot(tm_10x, features = mousegene) + labs(x = '', y = '')  + coord_flip() + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5))
#ggsave('mouse_tm10x.pdf', units = 'cm', width = 14, height = 10)
VlnPlot(tm_10x, features = mousegene)

DotPlot(tm_ss2, features = mousegene) + labs(x = '', y = '')  + coord_flip() + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5))
#ggsave('mouse_tmss2.pdf', units = 'cm', width = 14, height = 10)
VlnPlot(tm_ss2, features = 'Meis1')


#mouse <- merge(x = mca, y = tm_merged)
#head(mouse@meta.data)
#unique(mouse@meta.data$tissue)
#mouse <- ScaleData(object = mouse, vars.to.regress = c('nCount_RNA', 'orig.ident'))

#unique(mouse@meta.data$celltype_refined)
#mouse@meta.data$celltype_refined <- factor(mouse@meta.data$celltype_refined,
#                                           levels = c("Progenitors", "Granulocytopoietic cells", "Erythroids", "Neutrophils", "Monocytes", 
#                                                      "Macrophage", "DC", "pDC", "Basophil", "T cells", "NK cells", "B cells", "Plasma cells"))
#mouse@meta.data$celltype_refined <- mapvalues(mouse@meta.data$celltype_refined,
#                                              from = levels(mouse@meta.data$celltype_refined),
#                                              to = c("Progenitors", "Gran. prog.", "Erythroids", "Neutrophils", "Monocytes", 
#                                                     "Macrophages", "DC", "pDC", "Basophils", "T cells", "NK cells", "B cells", "Plasma cells"))
#Idents(mouse) <- 'celltype_refined'
#saveRDS(mouse, '../../mouse/mouse_merged.Rds')

mouse <- readRDS('../../mouse/mouse_merged.Rds')
DotPlot(mouse, features = mousegene) + 
  labs(x = '', y = '', title = '') + coord_flip() + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5))
#ggsave('mouse_merged.pdf', units = 'cm', width = 12, height = 10)



dotplot1 <- DotPlot(subset(mouse, tissue == 'BoneMarrow'), features = mousegene) + 
  labs(x = '', y = '', title = 'BoneMarrow')  + coord_flip() + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5))
dotplot2 <- DotPlot(subset(mouse, tissue == 'PeripheralBlood'), features = mousegene) + 
  labs(x = '', y = '', title = 'PeripheralBlood')  + coord_flip() + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5))
CombinePlots(plots = list(dotplot1, dotplot2), ncol = 1)
#ggsave('mouse_merged.byOrigin.pdf', units = 'cm', width = 12, height = 18)
dev.off()



head(hcl@meta.data)
DotPlot(hcl, features = humangene) + labs(x = '', y = '')  + coord_flip() + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5))
#ggsave('human_hcl.pdf', units = 'cm', width = 14, height = 10)
DotPlot(subset(hcl, sample == 'AdultBoneMarrow'), features = humangene) + labs(x = '', y = '')  + coord_flip() + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5))
DotPlot(subset(hcl, sample == 'AdultPeripheralBlood'), features = humangene) + labs(x = '', y = '')  + coord_flip() + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5))



### Server ###
hca <- readRDS('tmp/integrated.hca.alignment.filtered.Rds')
humangene <- c('HMGB3', 'SUMO3', 'SET', 'PCNA', 'KPNB1', 'MEIS1')
DotPlot(hca, features = humangene) + labs(x = '', y = '')  + coord_flip() + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5))
#ggsave('human_hca.pdf', units = 'cm', width = 14, height = 10)


hca@meta.data$sample <- 'AdultBoneMarrow'
human <- merge(x = hca, y = hcl)
head(human@meta.data)
human <- ScaleData(object = human, vars.to.regress = c('nCount_RNA', 'orig.ident'))

unique(human@meta.data$celltype_refined)


human@meta.data$celltype_refined <- mapvalues(human@meta.data$celltype_refined,
                                              from = unique(human@meta.data$celltype_refined),
                                              to = c("Monocytes", "pre-B cells", "T cells", "NK cells", "DC", "Progenitors", "Erythroids", "B cells", "pro-B cells", 
                                                     "CD8 T cells", "pDC", "Plasma cells", "Granulocyte progenitor", "pre-PC", "Platelet", "Plasma cells", "Neutrophils", "Macrophage"))

human@meta.data$celltype_refined <- factor(human@meta.data$celltype_refined,
                                           levels = c("Progenitors", "Granulocyte progenitor", "Erythroids", "Neutrophils", "Monocytes", "Macrophage", "DC", "pDC",
                                                      "T cells", "CD8 T cells", "NK cells", "pre-PC", "pro-B cells", "pre-B cells", "B cells", "Plasma cells", "Platelet"))
Idents(human) <- 'celltype_refined'
saveRDS(human, '../human_merged.Rds')


humangene <- c('HMGB3', 'SUMO3', 'SET', 'PCNA', 'KPNB1', 'MEIS1')
dotplot1 <- DotPlot(human, features = humangene) + 
  labs(x = '', y = '', title = '')  + coord_flip() + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5))
ggsave('../human_merged.ph1.pdf', units = 'cm', width = 14, height = 10)

dotplot1 <- DotPlot(subset(human, sample == 'AdultBoneMarrow'), features = humangene) + 
  labs(x = '', y = '', title = 'Bone marrow')  + coord_flip() + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5))
dotplot2 <- DotPlot(subset(human, sample == 'AdultPeripheralBlood'), features = humangene) + 
  labs(x = '', y = '', title = 'Peripheral blood')  + coord_flip() + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5))
CombinePlots(plots = list(dotplot1, dotplot2), ncol = 1)
ggsave('../human_merged.ph1.byOrigin.pdf', units = 'cm', width = 14, height = 20)



