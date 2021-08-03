library(Seurat)
library(ggplot2)

fishobj <- readRDS('../../../Model_species/zebrafish/Tang_Q_et_al/InDrops/tmp/zebrafish.harmony.bloodcells.Rds')
head(fishobj@meta.data); nrow(fishobj@meta.data) # 3301
data <- as.matrix(GetAssayData(fishobj, slot = 'data')); dim(data) # 22842  3301

conserved <- read.delim('../../2.genelists/fish/diopt_fishtomouse.usegenes.txt')
conserved[1:6, 1:10]; nrow(conserved) # 8713


library(Rtsne)
library(plyr)
library(RColorBrewer)

use_cols <- colorRampPalette(brewer.pal(8, "Accent"))(18)
use_cols <- c('HSCs' = use_cols[c(1)],
              'Erythroid' = use_cols[c(3)],
              'Neutrophil' = use_cols[c(4)],
              'Macrophage' = use_cols[c(6)],
              'NK/T cell' = use_cols[c(10)],
              'B cell' = use_cols[c(16)])


for (seed in c(210730021:210730025)){
  set.seed(seed)
  
  ### random sampling ###
  rand_bc <- rownames(fishobj@meta.data) # use all cells
  #length(rand_bc) 
  
  data_mini <- t(data[intersect(as.character(conserved$ZebrafishGeneID), rownames(data)), rand_bc])
  dim(data_mini) # 3301 8713
  label_mini <- fishobj@meta.data[rownames(data_mini), ]
  
  ### get tSNE coordinates ###
  tsne <- Rtsne(data_mini, dims = 3, pca_scale = F)
  scores <- as.data.frame(tsne$Y)
  rownames(scores) <- rownames(data_mini)
  colnames(scores) <- c('tSNE1', 'tSNE2', 'tSNE3')
  scores$celltype <- label_mini$celltype
  

  plt <- ggplot(scores, aes(x = tSNE1, y = tSNE2, colour = celltype)) + 
    geom_point() + 
    scale_colour_manual(name="Cell type", values = use_cols) + 
    theme_bw() + 
    theme(axis.text = element_text(colour = 'black'),
          axis.ticks = element_line(colour = 'black'),
          panel.grid = element_blank()) +
    labs(title =  ''); plt
  ggsave(paste0(c('./', seed, '.1-2.pdf'), collapse = ''), units = 'cm', width = 11, height = 8)
  
  plt <- ggplot(scores, aes(x = tSNE1, y = tSNE3, colour = celltype)) + 
    geom_point() + 
    scale_colour_manual(name="Cell type", values = use_cols) + 
    theme_bw() + 
    theme(axis.text = element_text(colour = 'black'),
          axis.ticks = element_line(colour = 'black'),
          panel.grid = element_blank()) +
    labs(title =  ''); plt
  ggsave(paste0(c('./', seed, '.1-3.pdf'), collapse = ''), units = 'cm', width = 11, height = 8)
  
  plt <- ggplot(scores, aes(x = tSNE2, y = tSNE3, colour = celltype)) + 
    geom_point() + 
    scale_colour_manual(name="Cell type", values = use_cols) + 
    theme_bw() + 
    theme(axis.text = element_text(colour = 'black'),
          axis.ticks = element_line(colour = 'black'),
          panel.grid = element_blank()) +
    labs(title =  ''); plt
  ggsave(paste0(c('./', seed, '.2-3.pdf'), collapse = ''), units = 'cm', width = 11, height = 8)
  
  write.table(data.frame(Barcode = rownames(scores), scores, check.names = F, check.rows = F), 
              paste0(c('./', seed, '.txt'), collapse = ''), quote = F, sep = '\t', row.names = F, col.names = T)
}


