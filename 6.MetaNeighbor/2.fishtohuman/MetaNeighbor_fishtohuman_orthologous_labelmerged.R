library(MetaNeighbor)
library(SummarizedExperiment)
#browseVignettes("SummarizedExperiment")

### Orthologous genes ###
ortho <- read.delim('../../2.genelists/fish/diopt_fishtohuman.usegenes.txt')
ortho$PredictionDerivedFrom <- NULL
head(ortho); nrow(ortho)
length(unique(ortho$ZebrafishGeneID)); length(unique(ortho$HumanSymbol)) # 12310

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


hca_pseudocell_flt <- data.frame(matrix(ncol = ncol(hca_pseudocell), nrow = nrow(ortho))); dim(hca_pseudocell_flt) # 6497 2425
colnames(hca_pseudocell_flt) <- colnames(hca_pseudocell)
rownames(hca_pseudocell_flt) <- as.character(ortho$HumanSymbol)
hca_pseudocell_flt[is.na(hca_pseudocell_flt)] <- 0
hca_pseudocell_flt[intersect(as.character(ortho$HumanSymbol), unique(rownames(hca_pseudocell), rownames(hca_pseudocell_flt))), ] <- 
  hca_pseudocell[intersect(as.character(ortho$HumanSymbol), unique(rownames(hca_pseudocell), rownames(hca_pseudocell_flt))), ]
hca_pseudocell_flt <- as.matrix(hca_pseudocell_flt)
hca_pseudocell_flt[1:4, 1:4]; dim(hca_pseudocell_flt) # 12310  2425
rownames(hca_pseudocell_flt) <- as.character(ortho$ZebrafishGeneID) # change "target" gene id to "ref" gene id
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
hcl_pseudocell_flt[1:4, 1:4]; dim(hcl_pseudocell_flt) # 12310  2153
rownames(hcl_pseudocell_flt) <- as.character(ortho$ZebrafishGeneID) # change "target" gene id to "ref" gene id
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
dim(fish_pseudocell_flt); fish_pseudocell_flt[1:4, 1:4] # 12310   327
dim(human_pseudocell_flt); human_pseudocell_flt[1:4, 1:4] # 12310  4578

identical(rownames(fish_pseudocell_flt), rownames(human_pseudocell_flt))
merged <- cbind(fish_pseudocell_flt, human_pseudocell_flt)
merged[1:4, 1:4]; dim(merged) # 12310  4905
merged_label <- rbind(fish_label, human_label)
identical(colnames(merged), rownames(merged_label))

merged_se <- SummarizedExperiment(assays = list(counts=merged), colData = merged_label); merged_se
remove(fish_pseudocell_flt, fish_label, human_pseudocell_flt, human_label, merged_label)



### Unsupervised 
#dir.create('celltype_merged')
var_genes <-  variableGenes(dat = merged_se, exp_labels = merged_se$study_id); length(var_genes) # merged - 883
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
pheatmap(celltype_NV[colnames(celltype_NV)[7:ncol(celltype_NV)], colnames(celltype_NV)[1:6]],
         border_color = NA, fontsize = 9, cluster_rows = F, cluster_cols = F, 
         breaks = c(0:20)*0.05, color = rev(cols(20)), 
         cellheight = 8, cellwidth = 8, treeheight_row = 5, treeheight_col = 5,
         filename = 'celltype_merged/MetaNeighbor_flytofish_orthologous.celltype_AUROC.oneway.pdf')

