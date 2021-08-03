library(MetaNeighbor)
library(SummarizedExperiment)
#browseVignettes("SummarizedExperiment")

### Orthologous genes ###
ortho <- read.delim('../../2.genelists/fish/diopt_fishtohuman.usegenes.txt')
ortho$PredictionDerivedFrom <- NULL
head(ortho); nrow(ortho)
length(unique(ortho$ZebrafishGeneID)); length(unique(ortho$HumanSymbol)) # 12044

### Reference data - fish ###
fish_pseudocell <- readRDS('../../../Model_species/zebrafish/Tang_Q_et_al/InDrops/tmp/celltype.pseudocell10.Rds')
fish_pseudocell[1:4, 1:4]; dim(fish_pseudocell) # 22842   327

fish_pseudocell_flt <- data.frame(matrix(ncol = ncol(fish_pseudocell), nrow = nrow(ortho))); dim(fish_pseudocell_flt) # 12044   327
colnames(fish_pseudocell_flt) <- colnames(fish_pseudocell)
rownames(fish_pseudocell_flt) <- as.character(ortho$ZebrafishGeneID)
fish_pseudocell_flt[is.na(fish_pseudocell_flt)] <- 0

length(intersect(as.character(ortho$ZebrafishGeneID), rownames(fish_pseudocell)))
fish_pseudocell_flt[intersect(as.character(ortho$ZebrafishGeneID), rownames(fish_pseudocell)), ] <- 
  fish_pseudocell[intersect(as.character(ortho$ZebrafishGeneID), rownames(fish_pseudocell)), ]
fish_pseudocell_flt <- as.matrix(fish_pseudocell_flt)
head(fish_pseudocell_flt); dim(fish_pseudocell_flt) # 12044   327

fish_label <- DataFrame(row.names = colnames(fish_pseudocell_flt), celltype = unlist(lapply(colnames(fish_pseudocell_flt), function (x) unlist(strsplit(as.character(x), split = ' pc'))[1] )) )
fish_label$study_id <- 'Fish'
head(fish_label); nrow(fish_label)


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


human_pseudocell_flt <- data.frame(matrix(ncol = ncol(human_pseudocell), nrow = nrow(ortho))); dim(human_pseudocell_flt) # 12044  2618 -> 2617
colnames(human_pseudocell_flt) <- colnames(human_pseudocell)
rownames(human_pseudocell_flt) <- as.character(ortho$HumanSymbol)
human_pseudocell_flt[is.na(human_pseudocell_flt)] <- 0

human_pseudocell_flt[intersect(as.character(ortho$HumanSymbol), rownames(human_pseudocell)), ] <- 
  human_pseudocell[intersect(as.character(ortho$HumanSymbol), rownames(human_pseudocell)), ]
human_pseudocell_flt <- as.matrix(human_pseudocell_flt)
human_pseudocell_flt[1:4, 1:4]; dim(human_pseudocell_flt) # 12044  2618 -> 2617

rownames(human_pseudocell_flt) <- as.character(ortho$ZebrafishGeneID) # change "target" gene id to "ref" gene id
human_pseudocell_flt[1:4, 1:4]

summary(as.factor(human_label_flt$celltype))
remove(human_pseudocell, human_label)


### 
dim(fish_pseudocell_flt); fish_pseudocell_flt[1:4, 1:4] # 12044   327
dim(human_pseudocell_flt); human_pseudocell_flt[1:4, 1:4] # 12044  2618 -> 2617

identical(rownames(fish_pseudocell_flt), rownames(human_pseudocell_flt))
merged <- cbind(fish_pseudocell_flt, human_pseudocell_flt)
merged[1:4, 1:4]; dim(merged) # 12044  2945
merged_label <- rbind(fish_label, human_label_flt)
identical(colnames(merged), rownames(merged_label))

merged_se <- SummarizedExperiment(assays = list(counts=merged), colData = merged_label); merged_se
remove(merged, fish_pseudocell_flt, human_pseudocell_flt, human_label_flt, merged_label)


### Unsupervised 
#dir.create('celltype_noPlt')
var_genes <-  variableGenes(dat = merged_se, exp_labels = merged_se$study_id); length(var_genes) # 881
#write.table(var_genes, 'celltype_noPlt/MetaNeighbor_fishtohuman_orthologous.celltype_vargenes.txt', quote = F, sep = '\n', row.names = F, col.names = F)

celltype_NV <- MetaNeighborUS(dat = merged_se,
                              study_id = merged_se$study_id,
                              cell_type = merged_se$celltype,
                              var_genes = var_genes)
celltype_NV_w <- celltype_NV[c( 1:length(unique(fish_label$celltype)) ), setdiff( c( 1:ncol(celltype_NV)), c(1:length(unique(fish_label$celltype)) ) )]
celltype_NV_w <- data.frame(celltype = rownames(celltype_NV_w), celltype_NV_w, check.names = F, check.rows = F)
#write.table(celltype_NV_w, 'celltype_noPlt/MetaNeighbor_fishtohuman_orthologous.celltype_NV.txt', quote = F, sep = '\t', row.names = F, col.names = T)

top_hits <- topHits(cell_NV = celltype_NV,
                    dat = merged_se,
                    study_id = merged_se$study_id,
                    cell_type = merged_se$celltype, 
                    threshold = 0.80); top_hits
#write.table(top_hits, 'celltype_noPlt/MetaNeighbor_fishtohuman_orthologous.celltype_top_hits.txt', quote = F, sep = '\t', row.names = F, col.names = T)

library(pheatmap)
library(RColorBrewer)

cols <- colorRampPalette(brewer.pal(11, "RdYlBu"))
pheatmap(celltype_NV, border_color = NA, fontsize = 9, cluster_rows = F, cluster_cols = F, 
         breaks = c(0:20)*0.05, color = rev(cols(20)), 
         cellheight = 8, cellwidth = 8, treeheight_row = 5, treeheight_col = 5,
         filename = 'celltype_noPlt/MetaNeighbor_fishtohuman_orthologous.celltype_AUROC.pdf')
dev.off()

colnames(celltype_NV)
pheatmap(celltype_NV[colnames(celltype_NV)[1:6], colnames(celltype_NV)[7:ncol(celltype_NV)]], 
         border_color = NA, fontsize = 9, cluster_rows = F, cluster_cols = F, 
         breaks = c(0:20)*0.05, color = rev(cols(20)), 
         cellheight = 8, cellwidth = 8, treeheight_row = 5, treeheight_col = 5,
         filename = 'celltype_noPlt/MetaNeighbor_fishtohuman_orthologous.celltype_AUROC.oneway.pdf')

