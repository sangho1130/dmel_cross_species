library(MetaNeighbor)
library(SummarizedExperiment)
#browseVignettes("SummarizedExperiment")

### Orthologous genes ###
ortho <- read.delim('../../2.genelists/fish/diopt_fishtomouse.usegenes.txt')
ortho$PredictionDerivedFrom <- NULL
head(ortho); nrow(ortho)
length(unique(ortho$ZebrafishGeneID)); length(unique(ortho$MouseSymbol)) # 12256

### Reference data - fish ###
fish_pseudocell <- readRDS('../../zebrafish/Tang_Q_et_al/InDrops/tmp/celltype.pseudocell10.Rds')
fish_pseudocell[1:4, 1:4]; dim(fish_pseudocell)

fish_pseudocell_flt <- data.frame(matrix(ncol = ncol(fish_pseudocell), nrow = nrow(ortho))); dim(fish_pseudocell_flt)
colnames(fish_pseudocell_flt) <- colnames(fish_pseudocell)
rownames(fish_pseudocell_flt) <- as.character(ortho$ZebrafishGeneID)
fish_pseudocell_flt[is.na(fish_pseudocell_flt)] <- 0

length(intersect(as.character(ortho$ZebrafishGeneID), rownames(fish_pseudocell)))
fish_pseudocell_flt[intersect(as.character(ortho$ZebrafishGeneID), rownames(fish_pseudocell)), ] <- 
  fish_pseudocell[intersect(as.character(ortho$ZebrafishGeneID), rownames(fish_pseudocell)), ]
fish_pseudocell_flt <- as.matrix(fish_pseudocell_flt)
head(fish_pseudocell_flt); dim(fish_pseudocell_flt)

fish_label <- DataFrame(row.names = colnames(fish_pseudocell_flt), celltype = unlist(lapply(colnames(fish_pseudocell_flt), function (x) unlist(strsplit(as.character(x), split = ' pc'))[1] )) )
fish_label$study_id <- 'Fish'
head(fish_label); nrow(fish_label)

### Target data - mouse, MCA ###
mouse_mca_pseudocell <- readRDS('../../mouse/mca/tmp/celltype_refined_noNeu.pseudocell10.Rds')
mouse_mca_pseudocell <- as.matrix(mouse_mca_pseudocell)
mouse_mca_pseudocell[1:4, 1:4]; dim(mouse_mca_pseudocell) # 15258  911

mouse_mca_pseudocell_flt <- data.frame(matrix(ncol = ncol(mouse_mca_pseudocell), nrow = nrow(ortho))); dim(mouse_mca_pseudocell_flt)
colnames(mouse_mca_pseudocell_flt) <- colnames(mouse_mca_pseudocell)
rownames(mouse_mca_pseudocell_flt) <- as.character(ortho$MouseSymbol)
mouse_mca_pseudocell_flt[is.na(mouse_mca_pseudocell_flt)] <- 0

mouse_mca_pseudocell_flt[intersect(as.character(ortho$MouseSymbol), rownames(mouse_mca_pseudocell)), ] <- mouse_mca_pseudocell[intersect(as.character(ortho$MouseSymbol), rownames(mouse_mca_pseudocell)), ]
mouse_mca_pseudocell_flt <- as.matrix(mouse_mca_pseudocell_flt)
head(mouse_mca_pseudocell_flt); dim(mouse_mca_pseudocell_flt)
rownames(mouse_mca_pseudocell_flt) <- as.character(ortho$ZebrafishGeneID) # change "target" gene id to "ref" gene id
colnames(mouse_mca_pseudocell_flt) <- unlist( lapply(colnames(mouse_mca_pseudocell_flt), function (x) gsub(as.character(x), pattern = 'pc', replacement = 'mca') ) )

mouse_mca_label <- DataFrame(row.names = colnames(mouse_mca_pseudocell_flt), 
                             celltype = unlist(lapply(colnames(mouse_mca_pseudocell_flt), function (x) unlist(strsplit(as.character(x), split = ' mca'))[1] )) )
mouse_mca_label$study_id <- 'Mouse'
head(mouse_mca_label); nrow(mouse_mca_label)


### Target data - mouse, Tabula Muris ###
mouse_tm_pseudocell <- readRDS('../../mouse/tabulamuris_10x/tmp/celltype_refined.pseudocell10.Rds')
mouse_tm_pseudocell <- as.matrix(mouse_tm_pseudocell)
mouse_tm_pseudocell[1:4, 1:4]; dim(mouse_tm_pseudocell)

mouse_tm_pseudocell_flt <- data.frame(matrix(ncol = ncol(mouse_tm_pseudocell), nrow = nrow(ortho))); dim(mouse_tm_pseudocell_flt)
colnames(mouse_tm_pseudocell_flt) <- colnames(mouse_tm_pseudocell)
rownames(mouse_tm_pseudocell_flt) <- as.character(ortho$MouseSymbol)
mouse_tm_pseudocell_flt[is.na(mouse_tm_pseudocell_flt)] <- 0

mouse_tm_pseudocell_flt[intersect(as.character(ortho$MouseSymbol), rownames(mouse_tm_pseudocell)), ] <- mouse_tm_pseudocell[intersect(as.character(ortho$MouseSymbol), rownames(mouse_tm_pseudocell)), ]
mouse_tm_pseudocell_flt <- as.matrix(mouse_tm_pseudocell_flt)
rownames(mouse_tm_pseudocell_flt) <- as.character(ortho$ZebrafishGeneID) # change "target" gene id to "ref" gene id
colnames(mouse_tm_pseudocell_flt) <- unlist( lapply(colnames(mouse_tm_pseudocell_flt), function (x) gsub(as.character(x), pattern = 'pc', replacement = 'tm10x') ) )

mouse_tm_label <- DataFrame(row.names = colnames(mouse_tm_pseudocell_flt), 
                            celltype = unlist(lapply(colnames(mouse_tm_pseudocell_flt), function (x) unlist(strsplit(as.character(x), split = ' tm10x'))[1] )) )
mouse_tm_label$study_id <- 'Mouse'
head(mouse_tm_label); nrow(mouse_tm_label)


### Target data - mouse, Tabula Muris SS2 ###
mouse_tmss2_pseudocell <- readRDS('../../mouse/tabulamuris_ss2/tmp/celltype_refined.pseudocell10.Rds')
mouse_tmss2_pseudocell <- as.matrix(mouse_tmss2_pseudocell)
mouse_tmss2_pseudocell[1:4, 1:4]; dim(mouse_tmss2_pseudocell)

mouse_tmss2_pseudocell_flt <- data.frame(matrix(ncol = ncol(mouse_tmss2_pseudocell), nrow = nrow(ortho))); dim(mouse_tmss2_pseudocell_flt)
colnames(mouse_tmss2_pseudocell_flt) <- colnames(mouse_tmss2_pseudocell)
rownames(mouse_tmss2_pseudocell_flt) <- as.character(ortho$MouseSymbol)
mouse_tmss2_pseudocell_flt[is.na(mouse_tmss2_pseudocell_flt)] <- 0

mouse_tmss2_pseudocell_flt[intersect(as.character(ortho$MouseSymbol), rownames(mouse_tmss2_pseudocell)), ] <- 
  mouse_tmss2_pseudocell[intersect(as.character(ortho$MouseSymbol), rownames(mouse_tmss2_pseudocell)), ]
mouse_tmss2_pseudocell_flt <- as.matrix(mouse_tmss2_pseudocell_flt)
rownames(mouse_tmss2_pseudocell_flt) <- as.character(ortho$ZebrafishGeneID) # change "target" gene id to "ref" gene id
colnames(mouse_tmss2_pseudocell_flt) <- unlist( lapply(colnames(mouse_tmss2_pseudocell_flt), function (x) gsub(as.character(x), pattern = 'pc', replacement = 'tmss2') ) )

mouse_tmss2_label <- DataFrame(row.names = colnames(mouse_tmss2_pseudocell_flt), 
                            celltype = unlist(lapply(colnames(mouse_tmss2_pseudocell_flt), function (x) unlist(strsplit(as.character(x), split = ' tmss2'))[1] )) )
mouse_tmss2_label$study_id <- 'Mouse'
head(mouse_tmss2_label); nrow(mouse_tmss2_label)


### 
identical(rownames(fish_pseudocell_flt), rownames(mouse_mca_pseudocell_flt))
identical(rownames(fish_pseudocell_flt), rownames(mouse_tm_pseudocell_flt))
identical(rownames(fish_pseudocell_flt), rownames(mouse_tmss2_pseudocell_flt))

merged <- cbind(fish_pseudocell_flt, mouse_mca_pseudocell_flt, mouse_tm_pseudocell_flt, mouse_tmss2_pseudocell_flt); dim(merged)
merged_label <- rbind(fish_label, mouse_mca_label, mouse_tm_label, mouse_tmss2_label)
merged_se <- SummarizedExperiment(assays = list(counts=merged), colData = merged_label); merged_se


### Unsupervised 
#dir.create('celltype_merged')
var_genes <-  variableGenes(dat = merged_se, exp_labels = merged_se$study_id); length(var_genes) # merged - 898
#write.table(var_genes, 'celltype_merged/MetaNeighbor_flytofish_orthologous.celltype_vargenes.txt', quote = F, sep = '\n', row.names = F, col.names = F)

celltype_NV <- MetaNeighborUS(dat = merged_se,
                              study_id = merged_se$study_id,
                              cell_type = merged_se$celltype,
                              var_genes = var_genes)
celltype_NV_w <- celltype_NV[c( 1:length(unique(fish_label$celltype)) ), setdiff( c( 1:ncol(celltype_NV)), c(1:length(unique(fish_label$celltype)) ) )]
celltype_NV_w <- data.frame(celltype = rownames(celltype_NV_w), celltype_NV_w, check.names = F, check.rows = F)
#write.table(celltype_NV_w, 'celltype_merged/MetaNeighbor_flytofish_orthologous.celltype_NV.txt', quote = F, sep = '\t', row.names = F, col.names = T)

top_hits <- topHits(cell_NV = celltype_NV,
                    dat = merged_se,
                    study_id = merged_se$study_id,
                    cell_type = merged_se$celltype, 
                    threshold = 0.80); top_hits
#write.table(top_hits, 'celltype_merged/MetaNeighbor_flytofish_orthologous.celltype_top_hits.txt', quote = F, sep = '\t', row.names = F, col.names = T)

library(pheatmap)
library(RColorBrewer)

cols <- colorRampPalette(brewer.pal(11, "RdYlBu"))
pheatmap(celltype_NV, border_color = NA, fontsize = 9, cluster_rows = F, cluster_cols = F, 
         breaks = c(0:20)*0.05, color = rev(cols(20)), 
         cellheight = 8, cellwidth = 8, treeheight_row = 5, treeheight_col = 5)#,
         #filename = 'celltype_merged/MetaNeighbor_flytofish_orthologous.celltype_AUROC.pdf')
dev.off()

colnames(celltype_NV)
pheatmap(celltype_NV[colnames(celltype_NV)[7:19], colnames(celltype_NV)[1:6]],
         border_color = NA, fontsize = 9, cluster_rows = F, cluster_cols = F, 
         breaks = c(0:20)*0.05, color = rev(cols(20)), 
         cellheight = 8, cellwidth = 8, treeheight_row = 5, treeheight_col = 5,
         filename = 'celltype_merged/MetaNeighbor_flytofish_orthologous.celltype_AUROC.oneway.pdf')

