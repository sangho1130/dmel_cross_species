library(MetaNeighbor)
library(SummarizedExperiment)
library(plyr)
#browseVignettes("SummarizedExperiment")

### Orthologous genes ###
ortho <- read.delim('../../2.genelists/mouse/diopt_mousetohuman.usegenes.txt')
head(ortho); nrow(ortho)
length(unique(ortho$MouseGeneID)); length(unique(ortho$HumanSymbol)) # 15840

### Reference data - mouse ###
### Target data - mouse, MCA ###
mouse_mca_pseudocell <- readRDS('../../mouse/mca/tmp/celltype_refined_noNeu.pseudocell10.Rds')
mouse_mca_pseudocell <- as.matrix(mouse_mca_pseudocell)
mouse_mca_pseudocell[1:4, 1:4]; dim(mouse_mca_pseudocell) # 15258  911

mouse_mca_pseudocell_flt <- data.frame(matrix(ncol = ncol(mouse_mca_pseudocell), nrow = nrow(ortho))); dim(mouse_mca_pseudocell_flt)
colnames(mouse_mca_pseudocell_flt) <- colnames(mouse_mca_pseudocell)
rownames(mouse_mca_pseudocell_flt) <- as.character(ortho$MouseGeneID)
mouse_mca_pseudocell_flt[is.na(mouse_mca_pseudocell_flt)] <- 0

mouse_mca_pseudocell_flt[intersect(as.character(ortho$MouseGeneID), rownames(mouse_mca_pseudocell)), ] <- 
  mouse_mca_pseudocell[intersect(as.character(ortho$MouseGeneID), rownames(mouse_mca_pseudocell)), ]
mouse_mca_pseudocell_flt <- as.matrix(mouse_mca_pseudocell_flt)
head(mouse_mca_pseudocell_flt); dim(mouse_mca_pseudocell_flt)
rownames(mouse_mca_pseudocell_flt) <- as.character(ortho$MouseGeneID) # change "target" gene id to "ref" gene id
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
rownames(mouse_tm_pseudocell_flt) <- as.character(ortho$MouseGeneID)
mouse_tm_pseudocell_flt[is.na(mouse_tm_pseudocell_flt)] <- 0

mouse_tm_pseudocell_flt[intersect(as.character(ortho$MouseGeneID), rownames(mouse_tm_pseudocell)), ] <- 
  mouse_tm_pseudocell[intersect(as.character(ortho$MouseGeneID), rownames(mouse_tm_pseudocell)), ]
mouse_tm_pseudocell_flt <- as.matrix(mouse_tm_pseudocell_flt)
rownames(mouse_tm_pseudocell_flt) <- as.character(ortho$MouseGeneID) # change "target" gene id to "ref" gene id
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
rownames(mouse_tmss2_pseudocell_flt) <- as.character(ortho$MouseGeneID)
mouse_tmss2_pseudocell_flt[is.na(mouse_tmss2_pseudocell_flt)] <- 0

mouse_tmss2_pseudocell_flt[intersect(as.character(ortho$MouseGeneID), rownames(mouse_tmss2_pseudocell)), ] <- 
  mouse_tmss2_pseudocell[intersect(as.character(ortho$MouseGeneID), rownames(mouse_tmss2_pseudocell)), ]
mouse_tmss2_pseudocell_flt <- as.matrix(mouse_tmss2_pseudocell_flt)
rownames(mouse_tmss2_pseudocell_flt) <- as.character(ortho$MouseGeneID) # change "target" gene id to "ref" gene id
colnames(mouse_tmss2_pseudocell_flt) <- unlist( lapply(colnames(mouse_tmss2_pseudocell_flt), function (x) gsub(as.character(x), pattern = 'pc', replacement = 'tmss2') ) )

mouse_tmss2_label <- DataFrame(row.names = colnames(mouse_tmss2_pseudocell_flt), 
                               celltype = unlist(lapply(colnames(mouse_tmss2_pseudocell_flt), function (x) unlist(strsplit(as.character(x), split = ' tmss2'))[1] )) )
mouse_tmss2_label$study_id <- 'Mouse'
head(mouse_tmss2_label); nrow(mouse_tmss2_label)

mouse_mca_pseudocell_flt[1:4, 1:4]; dim(mouse_mca_pseudocell_flt); head(mouse_mca_label, n=2)
mouse_tm_pseudocell_flt[1:4, 1:4]; dim(mouse_tm_pseudocell_flt); head(mouse_tm_label, n=2)
mouse_tmss2_pseudocell_flt[1:4, 1:4]; dim(mouse_tmss2_pseudocell_flt); head(mouse_tmss2_label, n=2)

mouse_pseudocell <- cbind(mouse_mca_pseudocell_flt, mouse_tm_pseudocell_flt, mouse_tmss2_pseudocell_flt)
mouse_pseudocell[1:4 ,1:4]; dim(mouse_pseudocell) # 15840  1771
mouse_label <- rbind(mouse_mca_label, mouse_tm_label, mouse_tmss2_label)
head(mouse_label)


### Target data - human ###
# HCA
hca_pseudocell <- readRDS('../../human/hca/tmp/celltype_refined.pseudocell10.Rds')
hca_pseudocell <- as.matrix(hca_pseudocell)
hca_pseudocell[1:4, 1:4]; dim(hca_pseudocell) # 41532 24334

hca_label <- DataFrame(row.names = colnames(hca_pseudocell), 
                       celltype = unlist(lapply(colnames(hca_pseudocell), function (x) unlist(strsplit(as.character(x), split = ' pc'))[1] )) )
hca_label$study_id <- 'Human'
trunc(summary(as.factor(hca_label$celltype))/10)
head(hca_label)

set.seed(21053101)
hca_label_flt <- data.frame(matrix(nrow = 0, ncol = 2))
colnames(hca_label_flt) <- colnames(hca_label)
for (ct in unique(hca_label$celltype)) {
  bcs <- sample(rownames(subset(hca_label, celltype == ct)), size = trunc(nrow(subset(hca_label, celltype == ct))/10))
  hca_label_flt <- rbind(hca_label_flt, hca_label[bcs, ])
}
head(hca_label_flt); nrow(hca_label_flt) # 2425
hca_pseudocell <- hca_pseudocell[, rownames(hca_label_flt)]
dim(hca_pseudocell) # 41532  2425


# HCL
hcl_pseudocell <- readRDS('../../human/hcl/tmp/celltype_refined.pseudocell10.Rds')
hcl_pseudocell <- as.matrix(hcl_pseudocell)
hcl_pseudocell[1:4, 1:4]; dim(hcl_pseudocell) # 19546  2153


hca_pseudocell_flt <- data.frame(matrix(ncol = ncol(hca_pseudocell), nrow = nrow(ortho))); dim(hca_pseudocell_flt)
colnames(hca_pseudocell_flt) <- colnames(hca_pseudocell)
rownames(hca_pseudocell_flt) <- as.character(ortho$HumanSymbol)
hca_pseudocell_flt[is.na(hca_pseudocell_flt)] <- 0
hca_pseudocell_flt[intersect(as.character(ortho$HumanSymbol), unique(rownames(hca_pseudocell), rownames(hca_pseudocell_flt))), ] <- 
  hca_pseudocell[intersect(as.character(ortho$HumanSymbol), unique(rownames(hca_pseudocell), rownames(hca_pseudocell_flt))), ]
hca_pseudocell_flt <- as.matrix(hca_pseudocell_flt)
hca_pseudocell_flt[1:4, 1:4]; dim(hca_pseudocell_flt) # 15840  2425
rownames(hca_pseudocell_flt) <- as.character(ortho$MouseGeneID) # change "target" gene id to "ref" gene id
colnames(hca_pseudocell_flt) <- unlist( lapply(colnames(hca_pseudocell_flt), function (x) gsub(as.character(x), pattern = 'pc', replacement = 'hca') ) )
rownames(hca_label_flt) <- colnames(hca_pseudocell_flt)
hca_pseudocell_flt[1:4, 1:4]


hcl_pseudocell_flt <- data.frame(matrix(ncol = ncol(hcl_pseudocell), nrow = nrow(ortho))); dim(hcl_pseudocell_flt)
colnames(hcl_pseudocell_flt) <- colnames(hcl_pseudocell)
rownames(hcl_pseudocell_flt) <- as.character(ortho$HumanSymbol)
hcl_pseudocell_flt[is.na(hcl_pseudocell_flt)] <- 0
hcl_pseudocell_flt[intersect(as.character(ortho$HumanSymbol), unique(rownames(hcl_pseudocell), rownames(hca_pseudocell_flt))), ] <- 
  hcl_pseudocell[intersect(as.character(ortho$HumanSymbol), unique(rownames(hcl_pseudocell), rownames(hca_pseudocell_flt))), ]
hcl_pseudocell_flt <- as.matrix(hcl_pseudocell_flt)
hcl_pseudocell_flt[1:4, 1:4]; dim(hcl_pseudocell_flt) # 15840  2153
rownames(hcl_pseudocell_flt) <- as.character(ortho$MouseGeneID) # change "target" gene id to "ref" gene id
colnames(hcl_pseudocell_flt) <- unlist( lapply(colnames(hcl_pseudocell_flt), function (x) gsub(as.character(x), pattern = 'pc', replacement = 'hcl') ) )
hcl_pseudocell_flt[1:4, 1:4]

remove(hcl_pseudocell); remove(hca_pseudocell); remove(hca_label)

hcl_label <- DataFrame(row.names = colnames(hcl_pseudocell_flt), 
                       celltype = unlist(lapply(colnames(hcl_pseudocell_flt), function (x) unlist(strsplit(as.character(x), split = ' hcl'))[1] )) )
hcl_label$study_id <- 'Human'
head(hcl_label); nrow(hcl_label)

summary(as.factor(hca_label_flt$celltype))
summary(as.factor(hcl_label$celltype))


human_pseudocell_flt <- cbind(hca_pseudocell_flt, hcl_pseudocell_flt)
human_label <- rbind(hca_label_flt, hcl_label)
identical(rownames(human_label), colnames(human_pseudocell_flt))

human_label$celltype <- mapvalues(human_label$celltype, from = unique(human_label$celltype),
                                  to = c("Progenitors", "Granulocyte progenitor", "Erythroids", "Monocytes", "DC", 
                                         "pDC", "T cells", "CD8 T cells", "NK cells", "pre-PC", 
                                         "pro-B cells", "pre-B cells", "B cells", "Plasma cells", "Platelet", 
                                         "Neutrophils", "Macrophage", "Plasma cells"))
summary(as.factor(human_label$celltype))
remove(hca_pseudocell_flt, hcl_pseudocell_flt, hca_label_flt, hcl_label)


### 
dim(mouse_pseudocell); mouse_pseudocell[1:4, 1:4]
dim(human_pseudocell_flt); human_pseudocell_flt[1:4, 1:4]

identical(rownames(mouse_pseudocell), rownames(human_pseudocell_flt))
merged <- cbind(mouse_pseudocell, human_pseudocell_flt); dim(merged)
merged_label <- rbind(mouse_label, human_label)
merged_se <- SummarizedExperiment(assays = list(counts=merged), colData = merged_label); merged_se
remove(mouse_mca_pseudocell, mouse_mca_pseudocell_flt, mouse_mca_label, mouse_pseudocell, 
       mouse_tm_pseudocell, mouse_tm_pseudocell_flt, mouse_tm_label, 
       mouse_tmss2_pseudocell, mouse_tmss2_pseudocell_flt, mouse_tmss2_label,
       human_pseudocell_flt, mouse_label, human_label, merged_label)

### Unsupervised 
#dir.create('celltype_merged')
var_genes <-  variableGenes(dat = merged_se, exp_labels = merged_se$study_id); length(var_genes) # 1422
#write.table(var_genes, 'celltype_merged/MetaNeighbor_orthologous.celltype_vargenes.txt', quote = F, sep = '\n', row.names = F, col.names = F)

celltype_NV <- MetaNeighborUS(dat = merged_se,
                              study_id = merged_se$study_id,
                              cell_type = merged_se$celltype,
                              var_genes = var_genes)
colnames(celltype_NV)
celltype_NV_w <- celltype_NV[c( 1:13 ), setdiff( c( 1:ncol(celltype_NV)), c(1:13 ) )]
celltype_NV_w <- data.frame(celltype = rownames(celltype_NV_w), celltype_NV_w, check.names = F, check.rows = F)
#write.table(celltype_NV_w, 'celltype_merged/MetaNeighbor_orthologous.celltype_NV.txt', quote = F, sep = '\t', row.names = F, col.names = T)

top_hits <- topHits(cell_NV = celltype_NV,
                    dat = merged_se,
                    study_id = merged_se$study_id,
                    cell_type = merged_se$celltype, 
                    threshold = 0.80); top_hits
#write.table(top_hits, 'celltype_merged/MetaNeighbor_orthologous.celltype_top_hits.txt', quote = F, sep = '\t', row.names = F, col.names = T)

library(pheatmap)
library(RColorBrewer)

cols <- colorRampPalette(brewer.pal(11, "RdYlBu"))
pheatmap(celltype_NV, border_color = NA, fontsize = 9,
         breaks = c(0:20)*0.05, color = rev(cols(20)), 
         cellheight = 8, cellwidth = 8, treeheight_row = 5, treeheight_col = 5)#,
         #filename = 'celltype_merged/MetaNeighbor_orthologous.celltype_AUROC.pdf')
dev.off()

colnames(celltype_NV)
pheatmap(celltype_NV[setdiff( c( 1:ncol(celltype_NV)), c(1:13 ) ), c( 1:13 )], 
         border_color = NA, fontsize = 9, cluster_rows = F, cluster_cols = F,
         breaks = c(0:20)*0.05, color = rev(cols(20)), 
         cellheight = 8, cellwidth = 8, treeheight_row = 5, treeheight_col = 5,
         filename = 'celltype_merged/MetaNeighbor_orthologous.celltype_AUROC.oneway.pdf')


