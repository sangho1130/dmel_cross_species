library(Seurat)
library(ggplot2)

mouse_mca <- readRDS('../../../mouse/mca/tmp/mca.harmony.clustered.v2.Rds')
mouse_mca <- subset(mouse_mca, celltype_refined != 'Neutrophils')
mouse_mca@meta.data <- droplevels(mouse_mca@meta.data)
head(mouse_mca@meta.data); nrow(mouse_mca@meta.data) # 9165
data <- as.matrix(GetAssayData(mouse_mca, slot = 'data')); dim(data)

conserved <- read.delim('../../../2.genelists/mouse/diopt_mousetohuman.usegenes.txt')
conserved[1:6, 1:10]; nrow(conserved) # 15840


library(Rtsne)
library(plyr)
library(RColorBrewer)

for (seed in c(210605621:210605625)){
  set.seed(seed)
  
  ### random sampling ###
  rand_bc <- c()
  for (ct in levels(mouse_mca@meta.data$celltype_refined)) {
    tmp <- subset(mouse_mca, celltype_refined == ct)
    use_bc <- sample(rownames(tmp@meta.data), size = round(nrow(tmp@meta.data)/5))
    rand_bc <- append(rand_bc, use_bc)
  }
  #length(rand_bc)
  
  data_mini <- t(data[intersect(as.character(conserved$MouseSymbol), rownames(data)), rand_bc])
  label_mini <- mouse_mca@meta.data[rownames(data_mini), ]
  #head(label_mini)
  
  ### get tSNE coordinates ###
  tsne <- Rtsne(data_mini, dims = 3, pca_scale = F)
  scores <- as.data.frame(tsne$Y)
  rownames(scores) <- rownames(data_mini)
  colnames(scores) <- c('tSNE1', 'tSNE2', 'tSNE3')
  scores$celltype <- label_mini$celltype_refined
  
  getPalette <- colorRampPalette(brewer.pal(8, "Set2"))
  plt <- ggplot(scores, aes(x = tSNE1, y = tSNE2, colour = celltype)) + 
    geom_point() + 
    scale_colour_manual(name="Cell type", values = getPalette(length(levels(scores$celltype)))) + 
    theme_bw() + 
    theme(axis.text = element_text(colour = 'black'),
          axis.ticks = element_line(colour = 'black'),
          panel.grid = element_blank()) +
    labs(title =  ''); plt
  ggsave(paste0(c('./', seed, '.1-2.pdf'), collapse = ''), units = 'cm', width = 11, height = 8)
  
  plt <- ggplot(scores, aes(x = tSNE1, y = tSNE3, colour = celltype)) + 
    geom_point() + 
    scale_colour_manual(name="Cell type", values = getPalette(length(levels(scores$celltype)))) + 
    theme_bw() + 
    theme(axis.text = element_text(colour = 'black'),
          axis.ticks = element_line(colour = 'black'),
          panel.grid = element_blank()) +
    labs(title =  ''); plt
  ggsave(paste0(c('./', seed, '.1-3.pdf'), collapse = ''), units = 'cm', width = 11, height = 8)
  
  plt <- ggplot(scores, aes(x = tSNE2, y = tSNE3, colour = celltype)) + 
    geom_point() + 
    scale_colour_manual(name="Cell type", values = getPalette(length(levels(scores$celltype)))) + 
    theme_bw() + 
    theme(axis.text = element_text(colour = 'black'),
          axis.ticks = element_line(colour = 'black'),
          panel.grid = element_blank()) +
    labs(title =  ''); plt
  ggsave(paste0(c('./', seed, '.2-3.pdf'), collapse = ''), units = 'cm', width = 11, height = 8)
  
  write.table(data.frame(Barcode = rownames(scores), scores, check.names = F, check.rows = F), 
              paste0(c('./', seed, '.txt'), collapse = ''), quote = F, sep = '\t', row.names = F, col.names = T)
}


