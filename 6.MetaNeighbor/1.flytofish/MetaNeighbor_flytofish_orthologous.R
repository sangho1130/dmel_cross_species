library(MetaNeighbor)
library(SummarizedExperiment)
#browseVignettes("SummarizedExperiment")

### Orthologous genes ###
ortho <- read.delim('../../2.genelists/fly/diopt_flytofish.usegenes.txt')
head(ortho); nrow(ortho)
length(unique(ortho$FlyBaseID)); length(unique(ortho$ZebrafishSymbol)) # 5676

### Reference data - fly ###
ref_pseudocell <- readRDS('../../fly/tmp/celltype_v2_pseudocell10.Rds')
ref_pseudocell <- readRDS('../../fly/tmp/celltype_v2_origin_pseudocell10.Rds')

ref_pseudocell[1:4, 1:4]; dim(ref_pseudocell)
ref_pseudocell_flt <- as.matrix(ref_pseudocell[intersect(as.character(ortho$FlyGeneID), rownames(ref_pseudocell)), ])
identical(rownames(ref_pseudocell_flt), as.character(ortho$FlyGeneID))

ref_label <- DataFrame(row.names = colnames(ref_pseudocell_flt), celltype = unlist(lapply(colnames(ref_pseudocell), function (x) unlist(strsplit(as.character(x), split = ' pc'))[1] )) )
ref_label$study_id <- 'Fly'
head(ref_label); nrow(ref_label) # 4357


### Target data - fish ###
fish_pseudocell <- readRDS('../../zebrafish/Tang_Q_et_al/InDrops/tmp/celltype.pseudocell10.Rds')
fish_pseudocell[1:4, 1:4]; dim(fish_pseudocell)

fish_pseudocell_flt <- data.frame(matrix(ncol = ncol(fish_pseudocell), nrow = nrow(ortho))); dim(fish_pseudocell_flt)
colnames(fish_pseudocell_flt) <- colnames(fish_pseudocell)
rownames(fish_pseudocell_flt) <- as.character(ortho$ZebrafishSymbol)
fish_pseudocell_flt[is.na(fish_pseudocell_flt)] <- 0

fish_pseudocell_flt[intersect(as.character(ortho$ZebrafishSymbol), rownames(fish_pseudocell)), ] <- fish_pseudocell[intersect(as.character(ortho$ZebrafishSymbol), rownames(fish_pseudocell)), ]
fish_pseudocell_flt <- as.matrix(fish_pseudocell_flt)
rownames(fish_pseudocell_flt) <- as.character(ortho$FlyGeneID) # change "target" gene id to "ref" gene id

fish_label <- DataFrame(row.names = colnames(fish_pseudocell_flt), celltype = unlist(lapply(colnames(fish_pseudocell_flt), function (x) unlist(strsplit(as.character(x), split = ' pc'))[1] )) )
fish_label$study_id <- 'Fish'
head(fish_label); nrow(fish_label)

merged <- cbind(ref_pseudocell_flt, fish_pseudocell_flt); dim(merged)
merged_label <- rbind(ref_label, fish_label)
head(merged_label)
merged_se <- SummarizedExperiment(assays = list(counts=merged), colData = merged_label); merged_se


### Unsupervised 
#dir.create('celltype')
#dir.create('celltype_origin')
var_genes <-  variableGenes(dat = merged_se, exp_labels = merged_se$study_id); length(var_genes) # 343 | 351
#write.table(var_genes, 'celltype/MetaNeighbor_flytofish_orthologous.celltype_vargenes.txt', quote = F, sep = '\n', row.names = F, col.names = F)
#write.table(var_genes, 'celltype_origin/MetaNeighbor_flytofish_orthologous.celltype_origin_vargenes.txt', quote = F, sep = '\n', row.names = F, col.names = F)

celltype_NV <- MetaNeighborUS(dat = merged_se,
                              study_id = merged_se$study_id,
                              cell_type = merged_se$celltype,
                              var_genes = var_genes)
celltype_NV_w <- celltype_NV[c( 1:length(unique(ref_label$celltype)) ), setdiff( c( 1:ncol(celltype_NV)), c(1:length(unique(ref_label$celltype)) ) )]
celltype_NV_w <- data.frame(celltype = rownames(celltype_NV_w), celltype_NV_w, check.names = F, check.rows = F)
#write.table(celltype_NV_w, 'celltype/MetaNeighbor_flytofish_orthologous.celltype_NV.txt', quote = F, sep = '\t', row.names = F, col.names = T)
#write.table(celltype_NV_w, 'celltype_origin/MetaNeighbor_flytofish_orthologous.celltype_origin_NV.txt', quote = F, sep = '\t', row.names = F, col.names = T)

top_hits <- topHits(cell_NV = celltype_NV,
                    dat = merged_se,
                    study_id = merged_se$study_id,
                    cell_type = merged_se$celltype, 
                    threshold = 0.75); top_hits
#write.table(top_hits, 'celltype/MetaNeighbor_flytofish_orthologous.celltype_top_hits.txt', quote = F, sep = '\t', row.names = F, col.names = T)
#write.table(top_hits, 'celltype_origin/MetaNeighbor_flytofish_orthologous.celltype_origin_top_hits.txt', quote = F, sep = '\t', row.names = F, col.names = T)

library(pheatmap)
library(RColorBrewer)

cols <- colorRampPalette(brewer.pal(11, "RdYlBu"))
pheatmap(celltype_NV, border_color = NA, fontsize = 9, cluster_rows = F, cluster_cols = F, 
         breaks = c(0:20)*0.05, color = rev(cols(20)), 
         cellheight = 8, cellwidth = 8, treeheight_row = 5, treeheight_col = 5)#,
         #filename = 'celltype/MetaNeighbor_flytofish_orthologous.celltype_AUROC.pdf')
         #filename = 'celltype_origin/MetaNeighbor_flytofish_orthologous.celltype_origin_AUROC.pdf')
dev.off()

pheatmap(celltype_NV, border_color = NA, fontsize = 9, 
         breaks = c(0:20)*0.05, color = rev(cols(20)), 
         cellheight = 8, cellwidth = 8, treeheight_row = 5, treeheight_col = 5)#,
         #filename = 'celltype/MetaNeighbor_flytofish_orthologous.celltype_AUROC.pdf')
         #filename = 'celltype_origin/MetaNeighbor_flytofish_orthologous.celltype_origin_AUROC.pdf')
dev.off()


#pheatmap(celltype_NV[ rownames(celltype_NV)[16:21], colnames(celltype_NV)[1:15]], 
pheatmap(celltype_NV[ rownames(celltype_NV)[9:14], colnames(celltype_NV)[1:8]], 
         border_color = NA, fontsize = 9, cluster_rows = F, cluster_cols = F, 
         breaks = c(0:20)*0.05, color = rev(cols(20)), 
         cellheight = 10, cellwidth = 10, treeheight_row = 5, treeheight_col = 5,
         filename = 'celltype/MetaNeighbor_flytofish_orthologous.celltype_AUROC.oneway.pdf')
         #filename = 'celltype_origin/MetaNeighbor_flytofish_orthologous.celltype_origin_AUROC.oneway.pdf')
