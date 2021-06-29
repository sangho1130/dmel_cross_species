library(Seurat)
library(ggplot2)

mouse_tm10x <- readRDS('../../../mouse/tabulamuris_10x/tmp/tm_10x.harmony.clustered.Rds')
head(mouse_tm10x@meta.data); nrow(mouse_tm10x@meta.data) # 3652
data <- as.matrix(GetAssayData(mouse_tm10x, slot = 'data')); dim(data)

conserved <- read.delim('../../../2.genelists/mouse/diopt_mousetohuman.usegenes.txt')
conserved[1:6, 1:10]; nrow(conserved) # 15840


library(Rtsne)
library(plyr)
library(RColorBrewer)

for (seed in c(210605721:210605725)){
  set.seed(seed)
  
  ### random sampling ###
  rand_bc <- c()
  for (ct in levels(mouse_tm10x@meta.data$celltype_refined)) {
    tmp <- subset(mouse_tm10x, celltype_refined == ct)
    use_bc <- sample(rownames(tmp@meta.data), size = round(nrow(tmp@meta.data)/5))
    rand_bc <- append(rand_bc, use_bc)
  }
  #length(rand_bc)
  
  data_mini <- t(data[intersect(as.character(conserved$MouseSymbol), rownames(data)), rand_bc])
  label_mini <- mouse_tm10x@meta.data[rownames(data_mini), ]
  
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
  ggsave(paste0(c('./', seed, '.1-2.pdf'), collapse = ''), units = 'cm', width = 12.5, height = 8)
  
  plt <- ggplot(scores, aes(x = tSNE1, y = tSNE3, colour = celltype)) + 
    geom_point() + 
    scale_colour_manual(name="Cell type", values = getPalette(length(levels(scores$celltype)))) + 
    theme_bw() + 
    theme(axis.text = element_text(colour = 'black'),
          axis.ticks = element_line(colour = 'black'),
          panel.grid = element_blank()) +
    labs(title =  ''); plt
  ggsave(paste0(c('./', seed, '.1-3.pdf'), collapse = ''), units = 'cm', width = 12.5, height = 8)
  
  plt <- ggplot(scores, aes(x = tSNE2, y = tSNE3, colour = celltype)) + 
    geom_point() + 
    scale_colour_manual(name="Cell type", values = getPalette(length(levels(scores$celltype)))) + 
    theme_bw() + 
    theme(axis.text = element_text(colour = 'black'),
          axis.ticks = element_line(colour = 'black'),
          panel.grid = element_blank()) +
    labs(title =  ''); plt
  ggsave(paste0(c('./', seed, '.2-3.pdf'), collapse = ''), units = 'cm', width = 12.5, height = 8)
  
  write.table(data.frame(Barcode = rownames(scores), scores, check.names = F, check.rows = F), 
              paste0(c('./', seed, '.txt'), collapse = ''), quote = F, sep = '\t', row.names = F, col.names = T)
}


