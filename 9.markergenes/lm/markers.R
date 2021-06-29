
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
flygene <- c('Arpc2', 'Arp2', 'Fim')
fishgene <- c('arpc2', 'actr2a', 'pls3')
mousegene <- c('Arpc2', 'Actr2', 'Lcp1')
humangene <- c('ARPC2', 'ACTR2', 'LCP1')


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
DotPlot(tm_ss2, features = mousegene) + labs(x = '', y = '')  + coord_flip() + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5))
#ggsave('mouse_tmss2.pdf', units = 'cm', width = 14, height = 10)


mouse <- readRDS('../../mouse/mouse_merged.Rds')
DotPlot(mouse, features = mousegene) + 
  labs(x = '', y = '', title = '')  + coord_flip() + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5))
#ggsave('mouse_merged.pdf', units = 'cm', width = 12, height = 9)

dotplot1 <- DotPlot(subset(mouse, tissue == 'BoneMarrow'), features = mousegene) + 
  labs(x = '', y = '', title = 'BoneMarrow')  + coord_flip() + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5))
dotplot2 <- DotPlot(subset(mouse, tissue == 'PeripheralBlood'), features = mousegene) + 
  labs(x = '', y = '', title = 'PeripheralBlood')  + coord_flip() + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5))
CombinePlots(plots = list(dotplot1, dotplot2), ncol = 1)
#ggsave('mouse_merged.byOrigin.pdf', units = 'cm', width = 12, height = 18)
dev.off()




DotPlot(hcl, features = humangene) + labs(x = '', y = '')  + coord_flip() + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5))
#ggsave('human_hcl.pdf', units = 'cm', width = 14, height = 10)



### Server ###
hca <- readRDS('tmp/integrated.hca.alignment.filtered.Rds')
humangene <- c('ARPC2', 'ACTR2', 'LCP1')
DotPlot(hca, features = humangene) + labs(x = '', y = '')  + coord_flip() + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5))
#ggsave('human_hca_lm.pdf', units = 'cm', width = 14, height = 10)


human <- readRDS('../human_merged.Rds')

humangene <- c('ARPC2', 'ACTR2', 'LCP1')
DotPlot(human, features = humangene) + 
  labs(x = '', y = '', title = '')  + coord_flip() + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5))
ggsave('../human_merged.lm.pdf', units = 'cm', width = 14, height = 10)

dotplot1 <- DotPlot(subset(human, sample == 'AdultBoneMarrow'), features = humangene) + 
  labs(x = '', y = '', title = 'Bone marrow')  + coord_flip() + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5))
dotplot2 <- DotPlot(subset(human, sample == 'AdultPeripheralBlood'), features = humangene) + 
  labs(x = '', y = '', title = 'Peripheral blood')  + coord_flip() + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5))
CombinePlots(plots = list(dotplot1, dotplot2), ncol = 1)
ggsave('../human_merged.lm.byOrigin.pdf', units = 'cm', width = 14, height = 20)

human_pb <- subset(human, sample == 'AdultPeripheralBlood')
head(human_pb@meta.data)
t.test(FetchData(human_pb, vars = 'LAPTM4A', cells = rownames(subset(human_pb@meta.data, celltype_refined %in% c('Monocytes', 'Macropahge'))) )[, 1] ,
       FetchData(human_pb, vars = 'LAPTM4A', cells = rownames(subset(human_pb@meta.data, !celltype_refined %in% c('Monocytes', 'Macropahge'))) )[, 1] ) # t = 17.31, df = 8948.6, p-value < 2.2e-16



