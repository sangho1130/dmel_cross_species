library(Seurat)
library(ggplot2)
library(Rtsne)
library(plyr)
library(RColorBrewer)

human <- readRDS('../../../Model_species/human/merged_commongenes/tmp/human.harmony_Library.flt.Rds')
human <- subset(human, celltype_refined != 'Platelet')
human@meta.data <- droplevels(human@meta.data)
head(human@meta.data); nrow(human@meta.data) # 262630


conserved <- read.delim('../../2.genelists/stats_orthologous.4219genes.txt')
head(conserved); nrow(conserved)


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



for (seed in c(210730821:210730830)){
  set.seed(seed)
  
  ### random sampling ###
  rand_bc <- c()
  for (ct in levels(human@meta.data$celltype_refined)) {
    tmp_meta <- subset(human@meta.data, celltype_refined == ct)
    use_bc <- sample(rownames(tmp_meta), size = round(nrow(tmp_meta)/50))
    rand_bc <- append(rand_bc, use_bc)
  }
  #length(rand_bc) # 4870
  
  data_mini <- as.matrix(FetchData(human, vars = as.character(conserved$HumanSymbol), cells = rand_bc, slot = 'data'))
  label_mini <- human@meta.data[rownames(data_mini), ]
  head(label_mini)
  
  ### get tSNE coordinates ###
  tsne <- Rtsne(data_mini, dims = 3, pca_scale = F)
  scores <- as.data.frame(tsne$Y)
  rownames(scores) <- rownames(data_mini)
  colnames(scores) <- c('tSNE1', 'tSNE2', 'tSNE3')
  scores$celltype_refined <- label_mini$celltype_refined
  
  getPalette <- colorRampPalette(brewer.pal(8, "Set2"))
  plt <- ggplot(scores, aes(x = tSNE1, y = tSNE2, colour = celltype_refined)) + 
    geom_point() + 
    scale_colour_manual(name="Cell type", values = use_cols) + 
    theme_bw() + 
    theme(axis.text = element_text(colour = 'black'),
          axis.ticks = element_line(colour = 'black'),
          panel.grid = element_blank()) +
    labs(title =  ''); plt
  ggsave(paste0(c('./', seed, '.1-2.pdf'), collapse = ''), units = 'cm', width = 12.5, height = 8)
  
  plt <- ggplot(scores, aes(x = tSNE1, y = tSNE3, colour = celltype_refined)) + 
    geom_point() + 
    scale_colour_manual(name="Cell type", values = use_cols) + 
    theme_bw() + 
    theme(axis.text = element_text(colour = 'black'),
          axis.ticks = element_line(colour = 'black'),
          panel.grid = element_blank()) +
    labs(title =  ''); plt
  ggsave(paste0(c('./', seed, '.1-3.pdf'), collapse = ''), units = 'cm', width = 12.5, height = 8)
  
  plt <- ggplot(scores, aes(x = tSNE2, y = tSNE3, colour = celltype_refined)) + 
    geom_point() + 
    scale_colour_manual(name="Cell type", values = use_cols) + 
    theme_bw() + 
    theme(axis.text = element_text(colour = 'black'),
          axis.ticks = element_line(colour = 'black'),
          panel.grid = element_blank()) +
    labs(title =  ''); plt
  ggsave(paste0(c('./', seed, '.2-3.pdf'), collapse = ''), units = 'cm', width = 12.5, height = 8)
  
  write.table(data.frame(Barcode = rownames(scores), scores, check.names = F, check.rows = F), 
              paste0(c('./', seed, '.txt'), collapse = ''), quote = F, sep = '\t', row.names = F, col.names = T)
}


