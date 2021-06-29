library(ggplot2)
library(plyr)
library(RColorBrewer)

tsne_coords <- list.files('./', pattern = '.txt')
celltypes <- unique(read.delim(tsne_coords[1], row.names = 1)$celltype)

use_cols <- colorRampPalette(brewer.pal(8, "Accent"))(18)
use_cols <- c('HSCs' = use_cols[c(1)],
              'Erythroid' = use_cols[c(3)],
              'Neutrophil' = use_cols[c(4)],
              'Macrophage' = use_cols[c(6)],
              'NK/T cell' = use_cols[c(10)],
              'B cell' = use_cols[c(16)])

for (afile in tsne_coords) {
  seed <- unlist(strsplit(afile, split = '\\.'))[1]
  tsne_coord <- read.delim(afile, row.names = 1)
  head(tsne_coord)
  
  plt <- ggplot(tsne_coord, aes(x = tSNE1, y = tSNE2, colour = celltype)) + 
    geom_point() + 
    scale_colour_manual(name="Cell type", values = use_cols) + 
    theme_bw() + 
    theme(axis.text = element_text(colour = 'black'),
          axis.ticks = element_line(colour = 'black'),
          panel.grid = element_blank()) +
    labs(title =  ''); plt
  ggsave(paste0(c('./', seed, '.1-2.pdf'), collapse = ''), units = 'cm', width = 11, height = 8)
  
  plt <- ggplot(tsne_coord, aes(x = tSNE1, y = tSNE3, colour = celltype)) + 
    geom_point() + 
    scale_colour_manual(name="Cell type", values = use_cols) + 
    theme_bw() + 
    theme(axis.text = element_text(colour = 'black'),
          axis.ticks = element_line(colour = 'black'),
          panel.grid = element_blank()) +
    labs(title =  ''); plt
  ggsave(paste0(c('./', seed, '.1-3.pdf'), collapse = ''), units = 'cm', width = 11, height = 8)
  
  plt <- ggplot(tsne_coord, aes(x = tSNE2, y = tSNE3, colour = celltype)) + 
    geom_point() + 
    scale_colour_manual(name="Cell type", values = use_cols) + 
    theme_bw() + 
    theme(axis.text = element_text(colour = 'black'),
          axis.ticks = element_line(colour = 'black'),
          panel.grid = element_blank()) +
    labs(title =  ''); plt
  ggsave(paste0(c('./', seed, '.2-3.pdf'), collapse = ''), units = 'cm', width = 11, height = 8)
  
}

