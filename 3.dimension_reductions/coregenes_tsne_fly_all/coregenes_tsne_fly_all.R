library(Seurat)
library(ggplot2)

flyobj <- readRDS('../../../Drop-seq_merge_harmony_6.22/tmp/dropseq.combined.harmony_flt.Rds')
head(flyobj@meta.data); nrow(flyobj@meta.data)
data <- as.matrix(GetAssayData(flyobj, slot = 'data')); dim(data)

conserved <- read.delim('../../2.genelists/stats_orthologous.4219genes.txt')
head(conserved)


library(Rtsne)
library(plyr)


for (seed in c(210730521:210730525)){
  set.seed(seed)
  
  ### random sampling ###
  rand_bc <- c()
  for (ct in levels(flyobj@meta.data$celltype)) {
    tmp <- subset(flyobj@meta.data, celltype == ct)
    use_bc <- sample(rownames(tmp), size = round(nrow(tmp)/10))
    rand_bc <- append(rand_bc, use_bc)
  }
  data_mini <- t(data[as.character(conserved$FlyGeneID), rand_bc])
  label_mini <- flyobj@meta.data[rownames(data_mini), ]
  
  ### get tSNE coordinates ###
  tsne <- Rtsne(data_mini, dims = 3, pca_scale = F)
  scores <- as.data.frame(tsne$Y)
  rownames(scores) <- rownames(data_mini)
  colnames(scores) <- c('tSNE1', 'tSNE2', 'tSNE3')
  scores$celltype <- label_mini$celltype
  
  plt <- ggplot(scores, aes(x = tSNE1, y = tSNE2, colour = celltype)) + 
    geom_point() + 
    scale_colour_manual(name="Cell type", 
                        values = c("PSC"="#f15fa6",
                                   "PH"="#207eb3",
                                   "PM"="#a80d0c",
                                   "LM"="#f0a142", 
                                   "CC"="#25a9b0", 
                                   "GST-rich" = "#a4a4a4",
                                   "Adipohemocyte"="#1a1a1a")) + 
    theme_bw() + 
    theme(axis.text = element_text(colour = 'black'),
          axis.ticks = element_line(colour = 'black'),
          panel.grid = element_blank()) +
    labs(title =  ''); plt
  ggsave(paste0(c('./', seed, '.1-2.pdf'), collapse = ''), units = 'cm', width = 11.5, height = 8)
  
  plt <- ggplot(scores, aes(x = tSNE1, y = tSNE3, colour = celltype)) + 
    geom_point() + 
    scale_colour_manual(name="Cell type", 
                        values = c("PSC"="#f15fa6",
                                   "PH"="#207eb3",
                                   "PM"="#a80d0c",
                                   "LM"="#f0a142", 
                                   "CC"="#25a9b0", 
                                   "GST-rich" = "#a4a4a4",
                                   "Adipohemocyte"="#1a1a1a")) + 
    theme_bw() + 
    theme(axis.text = element_text(colour = 'black'),
          axis.ticks = element_line(colour = 'black'),
          panel.grid = element_blank()) +
    labs(title =  ''); plt
  ggsave(paste0(c('./', seed, '.1-3.pdf'), collapse = ''), units = 'cm', width = 11.5, height = 8)
  
  plt <- ggplot(scores, aes(x = tSNE2, y = tSNE3, colour = celltype)) + 
    geom_point() + 
    scale_colour_manual(name="Cell type", 
                        values = c("PSC"="#f15fa6",
                                   "PH"="#207eb3",
                                   "PM"="#a80d0c",
                                   "LM"="#f0a142", 
                                   "CC"="#25a9b0", 
                                   "GST-rich" = "#a4a4a4",
                                   "Adipohemocyte"="#1a1a1a")) + 
    theme_bw() + 
    theme(axis.text = element_text(colour = 'black'),
          axis.ticks = element_line(colour = 'black'),
          panel.grid = element_blank()) +
    labs(title =  ''); plt
  ggsave(paste0(c('./', seed, '.2-3.pdf'), collapse = ''), units = 'cm', width = 11.5, height = 8)
  
  write.table(data.frame(Barcode = rownames(scores), scores, check.names = F, check.rows = F), 
              paste0(c('./', seed, '.txt'), collapse = ''), quote = F, sep = '\t', row.names = F, col.names = T)
}





