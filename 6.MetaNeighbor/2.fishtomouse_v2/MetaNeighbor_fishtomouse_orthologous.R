library(MetaNeighbor)
library(SummarizedExperiment)
#browseVignettes("SummarizedExperiment")

### Orthologous genes ###
ortho <- read.delim('../../2.genelists/fish/diopt_fishtomouse.usegenes.txt')
ortho$PredictionDerivedFrom <- NULL
head(ortho); nrow(ortho)
length(unique(ortho$ZebrafishGeneID)); length(unique(ortho$MouseSymbol)) # 8713

### Reference data - fish ###
fish_pseudocell <- readRDS('../../../Model_species/zebrafish/Tang_Q_et_al/InDrops/tmp/celltype.pseudocell10.Rds')
fish_pseudocell[1:4, 1:4]; dim(fish_pseudocell)

fish_pseudocell_flt <- data.frame(matrix(ncol = ncol(fish_pseudocell), nrow = nrow(ortho))); dim(fish_pseudocell_flt)
colnames(fish_pseudocell_flt) <- colnames(fish_pseudocell)
rownames(fish_pseudocell_flt) <- as.character(ortho$ZebrafishGeneID)
fish_pseudocell_flt[is.na(fish_pseudocell_flt)] <- 0

length(intersect(as.character(ortho$ZebrafishGeneID), rownames(fish_pseudocell)))
fish_pseudocell_flt[intersect(as.character(ortho$ZebrafishGeneID), rownames(fish_pseudocell)), ] <- 
  fish_pseudocell[intersect(as.character(ortho$ZebrafishGeneID), rownames(fish_pseudocell)), ]
fish_pseudocell_flt <- as.matrix(fish_pseudocell_flt)
head(fish_pseudocell_flt); dim(fish_pseudocell_flt) # 8713  327

fish_label <- DataFrame(row.names = colnames(fish_pseudocell_flt), celltype = unlist(lapply(colnames(fish_pseudocell_flt), function (x) unlist(strsplit(as.character(x), split = ' pc'))[1] )) )
fish_label$study_id <- 'Fish'
head(fish_label); nrow(fish_label)

### Target data - mouse ###
mouse_pseudocell <- readRDS('../../../Model_species/mouse/merged_commongenes_v3/tmp/celltype_refined.pseudocell10_v2.Rds')
mouse_pseudocell <- as.matrix(mouse_pseudocell)
mouse_pseudocell[1:4, 1:4]; dim(mouse_pseudocell) # 11504  1505

mouse_pseudocell_flt <- data.frame(matrix(ncol = ncol(mouse_pseudocell), nrow = nrow(ortho))); dim(mouse_pseudocell_flt)
colnames(mouse_pseudocell_flt) <- colnames(mouse_pseudocell)
rownames(mouse_pseudocell_flt) <- as.character(ortho$MouseSymbol)
mouse_pseudocell_flt[is.na(mouse_pseudocell_flt)] <- 0
dim(mouse_pseudocell_flt) # 8713 1505

length(intersect(as.character(ortho$MouseSymbol), rownames(mouse_pseudocell))) # 8713
mouse_pseudocell_flt[intersect(as.character(ortho$MouseSymbol), rownames(mouse_pseudocell)), ] <- 
  mouse_pseudocell[intersect(as.character(ortho$MouseSymbol), rownames(mouse_pseudocell)), ]
mouse_pseudocell_flt <- as.matrix(mouse_pseudocell_flt); dim(mouse_pseudocell_flt) # 8713 1505
rownames(mouse_pseudocell_flt) <- as.character(ortho$ZebrafishGeneID) # change "target" gene id to "ref" gene id
mouse_pseudocell_flt[1:4, 1:4]

mouse_label <- DataFrame(row.names = colnames(mouse_pseudocell_flt), 
                             celltype = unlist(lapply(colnames(mouse_pseudocell_flt), function (x) unlist(strsplit(as.character(x), split = ' pc'))[1] )) )
mouse_label$study_id <- 'Mouse'
head(mouse_label); nrow(mouse_label) # 1505


### 
identical(rownames(fish_pseudocell_flt), rownames(mouse_pseudocell_flt))

merged <- cbind(fish_pseudocell_flt, mouse_pseudocell_flt); dim(merged)
merged_label <- rbind(fish_label, mouse_label)
merged_se <- SummarizedExperiment(assays = list(counts=merged), colData = merged_label); merged_se


### Unsupervised 
#dir.create('celltype')
var_genes <-  variableGenes(dat = merged_se, exp_labels = merged_se$study_id); length(var_genes) # 551
#write.table(var_genes, 'celltype/MetaNeighbor_fishtomouse_orthologous.celltype_vargenes.txt', quote = F, sep = '\n', row.names = F, col.names = F)

celltype_NV <- MetaNeighborUS(dat = merged_se,
                              study_id = merged_se$study_id,
                              cell_type = merged_se$celltype,
                              var_genes = var_genes)
celltype_NV_w <- celltype_NV[c( 1:length(unique(fish_label$celltype)) ), setdiff( c( 1:ncol(celltype_NV)), c(1:length(unique(fish_label$celltype)) ) )]
celltype_NV_w <- data.frame(celltype = rownames(celltype_NV_w), celltype_NV_w, check.names = F, check.rows = F)
#write.table(celltype_NV_w, 'celltype/MetaNeighbor_fishtomouse_orthologous.celltype_NV.txt', quote = F, sep = '\t', row.names = F, col.names = T)

top_hits <- topHits(cell_NV = celltype_NV,
                    dat = merged_se,
                    study_id = merged_se$study_id,
                    cell_type = merged_se$celltype, 
                    threshold = 0.80); top_hits
#write.table(top_hits, 'celltype/MetaNeighbor_fishtomouse_orthologous.celltype_top_hits.txt', quote = F, sep = '\t', row.names = F, col.names = T)

library(pheatmap)
library(RColorBrewer)

cols <- colorRampPalette(brewer.pal(11, "RdYlBu"))
pheatmap(celltype_NV, border_color = NA, fontsize = 9, cluster_rows = F, cluster_cols = F, 
         breaks = c(0:20)*0.05, color = rev(cols(20)), 
         cellheight = 8, cellwidth = 8, treeheight_row = 5, treeheight_col = 5,
         filename = 'celltype/MetaNeighbor_fishtomouse_orthologous.celltype_AUROC.pdf')
#dev.off()

colnames(celltype_NV)
pheatmap(celltype_NV[colnames(celltype_NV)[1:6], colnames(celltype_NV)[7:19]],
         border_color = NA, fontsize = 9, cluster_rows = F, cluster_cols = F, 
         breaks = c(0:20)*0.05, color = rev(cols(20)), 
         cellheight = 8, cellwidth = 8, treeheight_row = 5, treeheight_col = 5,
         filename = 'celltype/MetaNeighbor_fishtomouse_orthologous.celltype_AUROC.oneway.pdf')

