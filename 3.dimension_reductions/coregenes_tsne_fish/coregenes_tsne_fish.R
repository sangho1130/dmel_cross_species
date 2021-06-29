library(Seurat)
library(ggplot2)

fishobj <- readRDS('../zebrafish/Tang_Q_et_al/InDrops/tmp/zebrafish.harmony.bloodcells.Rds')
head(fishobj@meta.data); nrow(fishobj@meta.data)
data <- as.matrix(GetAssayData(fishobj, slot = 'data')); dim(data)

conserved <- read.delim('../2.genelists/stats_orthologous.4847genes.txt')
head(conserved); nrow(conserved)


library(Rtsne)
library(plyr)
library(RColorBrewer)

dir.create('coregenes_tsne_fish')
for (seed in c(210524011:210524015)){
  set.seed(seed)
  
  ### random sampling ###
  #rand_bc <- c()
  #for (ct in levels(fishobj@meta.data$celltype)) {
  #  tmp <- subset(fishobj@meta.data, celltype == ct)
  #  use_bc <- sample(rownames(tmp), size = round(nrow(tmp)/10))
  #  rand_bc <- append(rand_bc, use_bc)
  #}
  rand_bc <- rownames(fishobj@meta.data) # use all cells
  data_mini <- t(data[as.character(conserved$ZebrafishGeneID), rand_bc]); dim(data_mini)
  label_mini <- fishobj@meta.data[rownames(data_mini), ]
  
  ### get tSNE coordinates ###
  tsne <- Rtsne(data_mini, dims = 3, pca_scale = F)
  scores <- as.data.frame(tsne$Y)
  rownames(scores) <- rownames(data_mini)
  colnames(scores) <- c('tSNE1', 'tSNE2', 'tSNE3')
  scores$celltype <- label_mini$celltype
  
  getPalette <- colorRampPalette(brewer.pal(8, "Set2"))
  plt <- ggplot(scores, aes(x = tSNE1, y = tSNE2, colour = celltype)) + 
    geom_point() + 
    scale_colour_manual(name="Cell type", values = getPalette(length(levels(fishobj@meta.data$celltype)))) + 
    theme_bw() + 
    theme(axis.text = element_text(colour = 'black'),
          axis.ticks = element_line(colour = 'black'),
          panel.grid = element_blank()) +
    labs(title =  ''); plt
  ggsave(paste0(c('coregenes_tsne_fish/', seed, '.1-2.pdf'), collapse = ''), units = 'cm', width = 11, height = 8)
  
  plt <- ggplot(scores, aes(x = tSNE1, y = tSNE3, colour = celltype)) + 
    geom_point() + 
    scale_colour_manual(name="Cell type", values = getPalette(length(levels(fishobj@meta.data$celltype)))) + 
    theme_bw() + 
    theme(axis.text = element_text(colour = 'black'),
          axis.ticks = element_line(colour = 'black'),
          panel.grid = element_blank()) +
    labs(title =  ''); plt
  ggsave(paste0(c('coregenes_tsne_fish/', seed, '.1-3.pdf'), collapse = ''), units = 'cm', width = 11, height = 8)
  
  plt <- ggplot(scores, aes(x = tSNE2, y = tSNE3, colour = celltype)) + 
    geom_point() + 
    scale_colour_manual(name="Cell type", values = getPalette(length(levels(fishobj@meta.data$celltype)))) + 
    theme_bw() + 
    theme(axis.text = element_text(colour = 'black'),
          axis.ticks = element_line(colour = 'black'),
          panel.grid = element_blank()) +
    labs(title =  ''); plt
  ggsave(paste0(c('coregenes_tsne_fish/', seed, '.2-3.pdf'), collapse = ''), units = 'cm', width = 11, height = 8)
  
  write.table(data.frame(Barcode = rownames(scores), scores, check.names = F, check.rows = F), 
              paste0(c('coregenes_tsne_fish/', seed, '.txt'), collapse = ''), quote = F, sep = '\t', row.names = F, col.names = T)
}



