library(Seurat)
library(ggplot2)

human_hca <- readRDS('../../human/hca/tmp/integrated.hca.alignment.filtered.Rds')
head(human_hca@meta.data); nrow(human_hca@meta.data) # 243398

#data <- as.matrix(GetAssayData(human_hca, slot = 'data')); dim(data)
#data <- readRDS('../../human/hca/tmp/expr.Rds')
#data <- log(data + 1); dim(data)

conserved <- read.delim('../../2.genelists/stats_orthologous.4847genes.txt')
head(conserved); nrow(conserved)


library(Rtsne)
library(plyr)
library(RColorBrewer)
for (seed in c(210615121:210615130)){
  set.seed(seed)
  
  ### random sampling ###
  rand_bc <- c()
  for (ct in levels(human_hca@meta.data$celltype_refined)) {
    tmp_meta <- subset(human_hca@meta.data, celltype_refined == ct)
    use_bc <- sample(rownames(tmp_meta), size = round(nrow(tmp_meta)/50))
    rand_bc <- append(rand_bc, use_bc)
  }
  #length(rand_bc) # 4870
  
  data_mini <- as.matrix(FetchData(human_hca, vars = as.character(conserved$HumanSymbol), cells = rand_bc, slot = 'data'))
  label_mini <- human_hca@meta.data[rownames(data_mini), ]
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
    scale_colour_manual(name="Cell type", values = getPalette(length(levels(human_hca@meta.data$celltype_refined)))) + 
    theme_bw() + 
    theme(axis.text = element_text(colour = 'black'),
          axis.ticks = element_line(colour = 'black'),
          panel.grid = element_blank()) +
    labs(title =  ''); plt
  ggsave(paste0(c('./', seed, '.1-2.pdf'), collapse = ''), units = 'cm', width = 12.5, height = 8)
  
  plt <- ggplot(scores, aes(x = tSNE1, y = tSNE3, colour = celltype_refined)) + 
    geom_point() + 
    scale_colour_manual(name="Cell type", values = getPalette(length(levels(human_hca@meta.data$celltype_refined)))) + 
    theme_bw() + 
    theme(axis.text = element_text(colour = 'black'),
          axis.ticks = element_line(colour = 'black'),
          panel.grid = element_blank()) +
    labs(title =  ''); plt
  ggsave(paste0(c('./', seed, '.1-3.pdf'), collapse = ''), units = 'cm', width = 12.5, height = 8)
  
  plt <- ggplot(scores, aes(x = tSNE2, y = tSNE3, colour = celltype_refined)) + 
    geom_point() + 
    scale_colour_manual(name="Cell type", values = getPalette(length(levels(human_hca@meta.data$celltype_refined)))) + 
    theme_bw() + 
    theme(axis.text = element_text(colour = 'black'),
          axis.ticks = element_line(colour = 'black'),
          panel.grid = element_blank()) +
    labs(title =  ''); plt
  ggsave(paste0(c('./', seed, '.2-3.pdf'), collapse = ''), units = 'cm', width = 12.5, height = 8)
  
  write.table(data.frame(Barcode = rownames(scores), scores, check.names = F, check.rows = F), 
              paste0(c('./', seed, '.txt'), collapse = ''), quote = F, sep = '\t', row.names = F, col.names = T)
}


