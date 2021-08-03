library(MetaNeighbor)
library(SummarizedExperiment)
library(plyr)
#browseVignettes("SummarizedExperiment")

### Orthologous genes ###
ortho <- read.delim('../../2.genelists/mouse/diopt_mousetohuman.usegenes.txt')
ortho$PredictionDerivedFrom <- NULL
head(ortho); nrow(ortho)
length(unique(ortho$MouseGeneID)); length(unique(ortho$HumanSymbol)) # 10380

### Reference data - mouse ###
mouse_pseudocell <- readRDS('../../../Model_species/mouse/merged_commongenes_v3/tmp/celltype_refined.pseudocell10_v2.Rds')
mouse_pseudocell <- as.matrix(mouse_pseudocell)
mouse_pseudocell[1:4, 1:4]; dim(mouse_pseudocell) # 11504  1505

mouse_pseudocell_flt <- data.frame(matrix(ncol = ncol(mouse_pseudocell), nrow = nrow(ortho))); dim(mouse_pseudocell_flt) # 10380  1505
colnames(mouse_pseudocell_flt) <- colnames(mouse_pseudocell)
rownames(mouse_pseudocell_flt) <- as.character(ortho$MouseGeneID)
mouse_pseudocell_flt[is.na(mouse_pseudocell_flt)] <- 0

mouse_pseudocell_flt[intersect(as.character(ortho$MouseGeneID), rownames(mouse_pseudocell)), ] <- 
  mouse_pseudocell[intersect(as.character(ortho$MouseGeneID), rownames(mouse_pseudocell)), ]
mouse_pseudocell_flt <- as.matrix(mouse_pseudocell_flt)
head(mouse_pseudocell_flt); dim(mouse_pseudocell_flt)
rownames(mouse_pseudocell_flt) <- as.character(ortho$MouseGeneID) # change "target" gene id to "ref" gene id

mouse_label <- DataFrame(row.names = colnames(mouse_pseudocell_flt), 
                             celltype = unlist(lapply(colnames(mouse_pseudocell_flt), function (x) unlist(strsplit(as.character(x), split = ' pc'))[1] )) )
mouse_label$study_id <- 'Mouse'
head(mouse_label); nrow(mouse_label)


### Target data - human ###
human_pseudocell <- readRDS('../../../Model_species/human/merged_commongenes/tmp/celltype_refined.pseudocell10.Rds')
human_pseudocell <- as.matrix(human_pseudocell)
human_pseudocell[1:4, 1:4]; dim(human_pseudocell) # 20131 26275

human_label <- DataFrame(row.names = colnames(human_pseudocell), 
                         celltype = unlist(lapply(colnames(human_pseudocell), function (x) unlist(strsplit(as.character(x), split = ' pc'))[1] )) )
human_label$study_id <- 'Human'
trunc(summary(as.factor(human_label$celltype))/10)
head(human_label)

set.seed(21072301)
human_label_flt <- data.frame(matrix(nrow = 0, ncol = 2))
colnames(human_label_flt) <- colnames(human_label)
for (ct in unique(human_label$celltype)[-17]) {
  bcs <- sample(rownames(subset(human_label, celltype == ct)), size = trunc(nrow(subset(human_label, celltype == ct))/10))
  human_label_flt <- rbind(human_label_flt, human_label[bcs, ])
}
head(human_label_flt); nrow(human_label_flt) # 2618 -> 2617
human_pseudocell <- human_pseudocell[, rownames(human_label_flt)]
human_pseudocell[1:4, 1:4]; dim(human_pseudocell) # 20131  2618 -> 2617


human_pseudocell_flt <- data.frame(matrix(ncol = ncol(human_pseudocell), nrow = nrow(ortho))); dim(human_pseudocell_flt) # 10380  2618 -> 2617
colnames(human_pseudocell_flt) <- colnames(human_pseudocell)
rownames(human_pseudocell_flt) <- as.character(ortho$HumanSymbol)
human_pseudocell_flt[is.na(human_pseudocell_flt)] <- 0

human_pseudocell_flt[intersect(as.character(ortho$HumanSymbol), rownames(human_pseudocell)), ] <- 
  human_pseudocell[intersect(as.character(ortho$HumanSymbol), rownames(human_pseudocell)), ]
human_pseudocell_flt <- as.matrix(human_pseudocell_flt)
human_pseudocell_flt[1:4, 1:4]; dim(human_pseudocell_flt) # 10380  2618 -> 2617

rownames(human_pseudocell_flt) <- as.character(ortho$MouseGeneID) # change "target" gene id to "ref" gene id
human_pseudocell_flt[1:4, 1:4]

summary(as.factor(human_label_flt$celltype))


### 
nrow(ortho)
dim(mouse_pseudocell_flt); mouse_pseudocell_flt[1:4, 1:4]
dim(human_pseudocell_flt); human_pseudocell_flt[1:4, 1:4]

identical(rownames(mouse_pseudocell_flt), rownames(human_pseudocell_flt))
merged <- cbind(mouse_pseudocell_flt, human_pseudocell_flt); dim(merged)
merged_label <- rbind(mouse_label, human_label_flt)
merged_se <- SummarizedExperiment(assays = list(counts=merged), colData = merged_label); merged_se
remove(mouse_pseudocell, mouse_pseudocell_flt, human_pseudocell, human_label, human_pseudocell_flt, merged)

### Unsupervised 
#dir.create('celltype_noPlt')
var_genes <-  variableGenes(dat = merged_se, exp_labels = merged_se$study_id); length(var_genes) # 769
write.table(var_genes, 'celltype_noPlt/MetaNeighbor_orthologous.celltype_vargenes.txt', quote = F, sep = '\n', row.names = F, col.names = F)

celltype_NV <- MetaNeighborUS(dat = merged_se,
                              study_id = merged_se$study_id,
                              cell_type = merged_se$celltype,
                              var_genes = var_genes)
celltype_NV_w <- celltype_NV[c( 1:13 ), c(14:nrow(celltype_NV)) ]
celltype_NV_w <- data.frame(celltype = rownames(celltype_NV_w), celltype_NV_w, check.names = F, check.rows = F)
#write.table(celltype_NV_w, 'celltype_noPlt/MetaNeighbor_orthologous.celltype_NV.txt', quote = F, sep = '\t', row.names = F, col.names = T)

top_hits <- topHits(cell_NV = celltype_NV,
                    dat = merged_se,
                    study_id = merged_se$study_id,
                    cell_type = merged_se$celltype,
                    threshold = 0.80); top_hits
#write.table(top_hits, 'celltype_noPlt/MetaNeighbor_orthologous.celltype_top_hits.txt', quote = F, sep = '\t', row.names = F, col.names = T)

library(pheatmap)
library(RColorBrewer)

cols <- colorRampPalette(brewer.pal(11, "RdYlBu"))
pheatmap(celltype_NV, border_color = NA, fontsize = 9,
         breaks = c(0:20)*0.05, color = rev(cols(20)), 
         cellheight = 8, cellwidth = 8, treeheight_row = 5, treeheight_col = 5,
         filename = 'celltype_noPlt/MetaNeighbor_orthologous.celltype_AUROC.pdf')
#dev.off()

colnames(celltype_NV)
pheatmap(celltype_NV[c( 1:13 ), c( 14:29 )], 
         border_color = NA, fontsize = 9, cluster_rows = F, cluster_cols = F,
         breaks = c(0:20)*0.05, color = rev(cols(20)), 
         cellheight = 8, cellwidth = 8, treeheight_row = 5, treeheight_col = 5,
         filename = 'celltype_noPlt/MetaNeighbor_orthologous.celltype_AUROC.oneway.pdf')


