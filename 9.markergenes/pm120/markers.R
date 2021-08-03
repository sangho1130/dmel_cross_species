
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


fish <- readRDS('../../../Model_species/zebrafish/Tang_Q_et_al/InDrops/tmp/zebrafish.harmony.bloodcells.Rds')
levels(Idents(fish))

mouse <- readRDS('../../../Model_species/mouse/merged_commongenes_v3/tmp/mouse.harmony.flt.Rds')
levels(Idents(mouse))

###
# /home/sangho/2020_newbuild/projects/drosophila/projectX_review/Model_species/human/merged_commongenes
human <- readRDS('../tmp/human.harmony_Library.flt.Rds')
human <- subset(human, celltype_refined != "Platelet")
human@meta.data <- droplevels(human@meta.data)
levels(Idents(human))

###
flygene <- c('Vps2', 'CG14767', 'Arpc3A', 'Spn55B', 'SkpA')
flygene <- c('Vps2', 'CG14767', 'Arpc3A')

fishgene <- c(NA, 'laptm4a', 'arpc3', 'serpinb1', 'skp1') 
fishgene <- c('laptm4a', 'arpc3')

mousegene <- c('Chmp2a', 'Laptm4a', 'Arpc3', 'Serpinb1a', 'Skp1a')
mousegene <- c('Chmp2a', 'Laptm4a', 'Arpc3')

humangene <- c('CHMP2A', 'LAPTM4A', 'ARPC3', 'SERPINB1', 'SKP1')
humangene <- c('CHMP2A', 'LAPTM4A', 'ARPC3')


DotPlot(bloodcells, features = flygene) + labs(x = '', y = '') + coord_flip() + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5))
#ggsave('drosophila_pm120.pdf', units = 'cm', width = 12, height = 10)


DotPlot(fish, features = fishgene) + labs(x = '', y = '')  + coord_flip() + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5))
#ggsave('zebrafish_pm120.pdf', units = 'cm', width = 11, height = 10)

fish_pseudoB <- readRDS('../../../Model_species/zebrafish/Tang_Q_et_al/InDrops/tmp/celltype.pseudo.Rds')
head(fish_pseudoB)
pheatmap::pheatmap(fish_pseudoB[fishgene, ], scale = 'row',
                   cellwidth = 10, cellheight = 10, cluster_rows = F, cluster_cols = F, filename = 'zebrafish_pm120.heatmap.pdf')


DotPlot(mouse, features = mousegene) + labs(x = '', y = '')  + coord_flip() + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5))
#ggsave('mouse_pm120.pdf', units = 'cm', width = 14, height = 11)

DotPlot(subset(mouse, tissue == 'BoneMarrow'), features = mousegene) + labs(x = '', y = '')  + 
  coord_flip() + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5))
#ggsave('mouse_pm120.bonemarrow.pdf', units = 'cm', width = 14, height = 11)
DotPlot(subset(mouse, tissue != 'BoneMarrow'), features = mousegene) + labs(x = '', y = '')  + 
  coord_flip() + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5))
#ggsave('mouse_pm120.peripheralblood.pdf', units = 'cm', width = 14, height = 11)

mouse_pseudoB <- readRDS('../../../Model_species/mouse/merged_commongenes_v3/tmp/mouse_pseudobulk_v2.Rds')
head(mouse_pseudoB)
pheatmap::pheatmap(mouse_pseudoB[mousegene, ], scale = 'row',
                   cellwidth = 10, cellheight = 10, cluster_rows = F, cluster_cols = F, filename = 'mouse_pm120.heatmap.pdf')

mouse_pseudoB_bm <- readRDS('../../../Model_species/mouse/merged_commongenes_v3/tmp/mouse_pseudobulk_v2_bonemarrow.Rds')
head(mouse_pseudoB_bm)
pheatmap::pheatmap(mouse_pseudoB_bm[mousegene, ], scale = 'row',
                   cellwidth = 10, cellheight = 10, cluster_rows = F, cluster_cols = F, filename = 'mouse_pm120.heatmap.bonemarrow.pdf')
mouse_pseudoB_pb <- readRDS('../../../Model_species/mouse/merged_commongenes_v3/tmp/mouse_pseudobulk_v2_peripheralblood.Rds')
head(mouse_pseudoB_pb)
pheatmap::pheatmap(mouse_pseudoB_pb[mousegene, ], scale = 'row',
                   cellwidth = 10, cellheight = 10, cluster_rows = F, cluster_cols = F, filename = 'mouse_pm120.heatmap.peripheralblood.pdf')



DotPlot(human, features = humangene) + labs(x = '', y = '') + coord_flip() + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5))
#ggsave('human_pm120.pdf', units = 'cm', width = 16, height = 11)

DotPlot(subset(human, sample == 'AdultBoneMarrow'), features = humangene) + labs(x = '', y = '') + 
  coord_flip() + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5))
ggsave('human_pm120.bonemarrow.pdf', units = 'cm', width = 16, height = 11)
DotPlot(subset(human, sample != 'AdultBoneMarrow'), features = humangene) + labs(x = '', y = '') + 
  coord_flip() + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5))
ggsave('human_pm120.peripheralblood.pdf', units = 'cm', width = 16, height = 11)


human_pseudoB <- readRDS('../../../Model_species/human/merged_commongenes/tmp/expr.pseudobulk.Rds')
head(human_pseudoB)
pheatmap::pheatmap(human_pseudoB[humangene, -17], scale = 'row',
                   cellwidth = 10, cellheight = 10, cluster_rows = F, cluster_cols = F, filename = 'human_pm120.heatmap.pdf')

human_pseudoB_bm <- readRDS('../../../Model_species/human/merged_commongenes/tmp/expr.pseudobulk_bonemarrow.Rds')
head(human_pseudoB_bm)
pheatmap::pheatmap(human_pseudoB_bm[humangene, -17], scale = 'row',
                   cellwidth = 10, cellheight = 10, cluster_rows = F, cluster_cols = F, filename = 'human_pm120.heatmap.bonemarrow.pdf')
human_pseudoB_pb <- readRDS('../../../Model_species/human/merged_commongenes/tmp/expr.pseudobulk_peripheralblood.Rds')
head(human_pseudoB_pb)
pheatmap::pheatmap(human_pseudoB_pb[humangene, -17], scale = 'row',
                   cellwidth = 10, cellheight = 10, cluster_rows = F, cluster_cols = F, filename = 'human_pm120.heatmap.peripheralblood.pdf')

