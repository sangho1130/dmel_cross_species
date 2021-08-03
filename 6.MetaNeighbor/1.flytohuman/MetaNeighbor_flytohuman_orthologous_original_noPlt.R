library(MetaNeighbor)
library(SummarizedExperiment)
#browseVignettes("SummarizedExperiment")

### Orthologous genes ###
ortho <- read.delim('../../2.genelists/fly/diopt_flytohuman.usegenes.txt')
head(ortho); nrow(ortho)
length(unique(ortho$FlyBaseID)); length(unique(ortho$HumanSymbol)) # 6385

### Reference data - fly ###
ref_pseudocell <- readRDS('../../../Model_species/fly/tmp/celltype_v2_pseudocell10.Rds')
ref_pseudocell <- readRDS('../../../Model_species/fly/tmp/celltype_v2_origin_pseudocell10.Rds')

ref_pseudocell[1:4, 1:4]; dim(ref_pseudocell)
ref_pseudocell_flt <- as.matrix(ref_pseudocell[intersect(as.character(ortho$FlyGeneID), rownames(ref_pseudocell)), ])
identical(rownames(ref_pseudocell_flt), as.character(ortho$FlyGeneID))

ref_label <- DataFrame(row.names = colnames(ref_pseudocell_flt), celltype = unlist(lapply(colnames(ref_pseudocell), function (x) unlist(strsplit(as.character(x), split = ' pc'))[1] )) )
ref_label$study_id <- 'Fly'
head(ref_label); nrow(ref_label) # 4354


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
dim(human_pseudocell) # 20131  2618 -> 2617


human_pseudocell_flt <- data.frame(matrix(ncol = ncol(human_pseudocell), nrow = nrow(ortho))); dim(human_pseudocell_flt) # 6385 2618 -> 2617
colnames(human_pseudocell_flt) <- colnames(human_pseudocell)
rownames(human_pseudocell_flt) <- as.character(ortho$HumanSymbol)
human_pseudocell_flt[is.na(human_pseudocell_flt)] <- 0

human_pseudocell_flt[intersect(as.character(ortho$HumanSymbol), rownames(human_pseudocell)), ] <- 
  human_pseudocell[intersect(as.character(ortho$HumanSymbol), rownames(human_pseudocell)), ]
human_pseudocell_flt <- as.matrix(human_pseudocell_flt)
human_pseudocell_flt[1:4, 1:4]; dim(human_pseudocell_flt) # 6385 2618 -> 2617
rownames(human_pseudocell_flt) <- as.character(ortho$FlyGeneID) # change "target" gene id to "ref" gene id
rownames(human_label_flt) <- colnames(human_pseudocell_flt)
human_pseudocell_flt[1:4, 1:4]

summary(as.factor(human_label$celltype))
remove(human_pseudocell, human_label)


### 
dim(ref_pseudocell_flt); ref_pseudocell_flt[1:4, 1:4]
dim(human_pseudocell_flt); human_pseudocell_flt[1:4, 1:4]

identical(rownames(ref_pseudocell_flt), rownames(human_pseudocell_flt))
merged <- cbind(ref_pseudocell_flt, human_pseudocell_flt)
merged[1:4, 1:4]; dim(merged) # 6385 6974 | 6971
merged_label <- rbind(ref_label, human_label_flt)
merged_se <- SummarizedExperiment(assays = list(counts=merged), colData = merged_label); merged_se
remove(ref_pseudocell, ref_pseudocell_flt, human_pseudocell_flt)



### Unsupervised 
#dir.create('celltype_noPlt')
#dir.create('celltype_origin_noPlt')
var_genes <-  variableGenes(dat = merged_se, exp_labels = merged_se$study_id); length(var_genes) # 344 | 345
#write.table(var_genes, 'celltype_noPlt/MetaNeighbor_flytohuman_orthologous.celltype_vargenes.txt', quote = F, sep = '\n', row.names = F, col.names = F)
#write.table(var_genes, 'celltype_origin_noPlt/MetaNeighbor_flytohuman_orthologous.celltype_origin_vargenes.txt', quote = F, sep = '\n', row.names = F, col.names = F)

celltype_NV <- MetaNeighborUS(dat = merged_se,
                              study_id = merged_se$study_id,
                              cell_type = merged_se$celltype,
                              var_genes = var_genes)
celltype_NV_w <- celltype_NV[c( 1:length(unique(ref_label$celltype)) ), setdiff( c( 1:ncol(celltype_NV)), c(1:length(unique(ref_label$celltype)) ) )]
celltype_NV_w <- data.frame(celltype = rownames(celltype_NV_w), celltype_NV_w, check.names = F, check.rows = F)
#write.table(celltype_NV_w, 'celltype_noPlt/MetaNeighbor_flytohuman_orthologous.celltype_NV.txt', quote = F, sep = '\t', row.names = F, col.names = T)
#write.table(celltype_NV_w, 'celltype_origin_noPlt/MetaNeighbor_flytohuman_orthologous.celltype_origin_NV.txt', quote = F, sep = '\t', row.names = F, col.names = T)

top_hits <- topHits(cell_NV = celltype_NV,
                    dat = merged_se,
                    study_id = merged_se$study_id,
                    cell_type = merged_se$celltype, 
                    threshold = 0.75); top_hits
#write.table(top_hits, 'celltype_noPlt/MetaNeighbor_flytohuman_orthologous.celltype_top_hits.txt', quote = F, sep = '\t', row.names = F, col.names = T)
#write.table(top_hits, 'celltype_origin_noPlt/MetaNeighbor_flytohuman_orthologous.celltype_origin_top_hits.txt', quote = F, sep = '\t', row.names = F, col.names = T)

library(pheatmap)
library(RColorBrewer)

cols <- colorRampPalette(brewer.pal(11, "RdYlBu"))
pheatmap(celltype_NV, border_color = NA, fontsize = 9, cluster_rows = F, cluster_cols = F, 
         breaks = c(0:20)*0.05, color = rev(cols(20)), 
         cellheight = 8, cellwidth = 8, treeheight_row = 5, treeheight_col = 5)#,
         #filename = 'celltype_noPlt/MetaNeighbor_flytohuman_orthologous.celltype_AUROC.pdf')
         #filename = 'celltype_origin_noPlt/MetaNeighbor_flytohuman_orthologous.celltype_origin_AUROC.pdf')
dev.off()

pheatmap(celltype_NV, border_color = NA, fontsize = 9, 
         breaks = c(0:20)*0.05, color = rev(cols(20)), 
         cellheight = 8, cellwidth = 8, treeheight_row = 5, treeheight_col = 5)#,
         #filename = 'celltype_noPlt/MetaNeighbor_flytohuman_orthologous.celltype_AUROC.clustered.pdf')
         #filename = 'celltype_origin_noPlt/MetaNeighbor_flytohuman_orthologous.celltype_origin_AUROC.clustered.pdf')
dev.off()


#pheatmap(celltype_NV[ rownames(celltype_NV)[16:nrow(celltype_NV)], colnames(celltype_NV)[1:15]], 
pheatmap(celltype_NV[ rownames(celltype_NV)[9:nrow(celltype_NV)], colnames(celltype_NV)[1:8]], 
         border_color = NA, fontsize = 9, cluster_rows = F, cluster_cols = F, 
         breaks = c(0:20)*0.05, color = rev(cols(20)), 
         cellheight = 10, cellwidth = 10, treeheight_row = 5, treeheight_col = 5)#,
         #filename = 'celltype/MetaNeighbor_flytohuman_orthologous.celltype_AUROC.oneway.pdf')
         #filename = 'celltype_origin_noPlt/MetaNeighbor_flytohuman_orthologous.celltype_origin_AUROC.oneway.pdf')
