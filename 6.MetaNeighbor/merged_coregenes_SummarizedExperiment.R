library(MetaNeighbor)
library(SummarizedExperiment)

### Orthologous genes ###
ortho <- read.delim('../2.genelists/stats_orthologous.4219genes.txt')
head(ortho)

### Reference data - fly ###
ref_pseudocell <- readRDS('../../Model_species/fly/tmp/celltype_v2_pseudocell10.Rds')
ref_pseudocell <- as.matrix(ref_pseudocell)
ref_pseudocell[1:4, 1:4]; dim(ref_pseudocell) # 16671  4357
ref_pseudocell_flt <- ref_pseudocell[intersect(as.character(ortho$FlyGeneID), rownames(ref_pseudocell)), ]
identical(rownames(ref_pseudocell_flt), as.character(ortho$FlyGeneID))

ref_label <- DataFrame(row.names = colnames(ref_pseudocell_flt), celltype = unlist(lapply(colnames(ref_pseudocell), function (x) unlist(strsplit(as.character(x), split = ' pc'))[1] )) )
ref_label$study_id <- 'Fly'
summary(as.factor(ref_label$celltype))
head(ref_label); nrow(ref_label) # 4357


### Target data - fish ###
fish_pseudocell <- readRDS('../../Model_species/zebrafish/Tang_Q_et_al/InDrops/tmp/celltype.pseudocell10.Rds')
fish_pseudocell <- as.matrix(fish_pseudocell)
fish_pseudocell[1:4, 1:4]; dim(fish_pseudocell)

fish_pseudocell_flt <- data.frame(matrix(ncol = ncol(fish_pseudocell), nrow = nrow(ortho))); dim(fish_pseudocell_flt)
colnames(fish_pseudocell_flt) <- colnames(fish_pseudocell)
rownames(fish_pseudocell_flt) <- as.character(ortho$ZebrafishGeneID)
fish_pseudocell_flt[is.na(fish_pseudocell_flt)] <- 0

fish_pseudocell_flt[intersect(as.character(ortho$ZebrafishGeneID), rownames(fish_pseudocell)), ] <- fish_pseudocell[intersect(as.character(ortho$ZebrafishGeneID), rownames(fish_pseudocell)), ]
fish_pseudocell_flt <- as.matrix(fish_pseudocell_flt)
rownames(fish_pseudocell_flt) <- as.character(ortho$FlyGeneID) # change "target" gene id to "ref" gene id

fish_label <- DataFrame(row.names = colnames(fish_pseudocell_flt), celltype = unlist(lapply(colnames(fish_pseudocell_flt), function (x) unlist(strsplit(as.character(x), split = ' pc'))[1] )) )
fish_label$study_id <- 'Fish'
head(fish_label); nrow(fish_label) # 327


### Target data - mouse ###
mouse_pseudocell <- readRDS('../../Model_species/mouse/merged_commongenes/tmp/celltype_refined.pseudocell10.Rds')
mouse_pseudocell[1:4, 1:4]; dim(mouse_pseudocell) # 11504  1672

mouse_pseudocell_flt <- data.frame(matrix(ncol = ncol(mouse_pseudocell), nrow = nrow(ortho))); dim(mouse_pseudocell_flt)
colnames(mouse_pseudocell_flt) <- colnames(mouse_pseudocell)
rownames(mouse_pseudocell_flt) <- as.character(ortho$MouseGeneID)
mouse_pseudocell_flt[is.na(mouse_pseudocell_flt)] <- 0

mouse_pseudocell_flt[intersect(as.character(ortho$MouseGeneID), rownames(mouse_pseudocell)), ] <- mouse_pseudocell[intersect(as.character(ortho$MouseGeneID), rownames(mouse_pseudocell)), ]
mouse_pseudocell_flt <- as.matrix(mouse_pseudocell_flt)
rownames(mouse_pseudocell_flt) <- as.character(ortho$FlyGeneID) # change "target" gene id to "ref" gene id

mouse_label <- DataFrame(row.names = colnames(mouse_pseudocell_flt), 
                             celltype = unlist(lapply(colnames(mouse_pseudocell_flt), function (x) unlist(strsplit(as.character(x), split = ' pc'))[1] )) )
mouse_label$study_id <- 'Mouse'
head(mouse_label); nrow(mouse_label) # 1672


### Target data - human  ###
human_pseudocell <- readRDS('../../Model_species/human/merged_commongenes/tmp/celltype_refined.pseudocell10.Rds')
human_pseudocell <- as.matrix(human_pseudocell)
human_pseudocell[1:4, 1:4]; dim(human_pseudocell) # 20131 26275

human_pseudocell_flt <- data.frame(matrix(ncol = ncol(human_pseudocell), nrow = nrow(ortho))); dim(human_pseudocell_flt)
colnames(human_pseudocell_flt) <- colnames(human_pseudocell)
rownames(human_pseudocell_flt) <- as.character(ortho$HumanSymbol)
human_pseudocell_flt[is.na(human_pseudocell_flt)] <- 0

human_pseudocell_flt[intersect(as.character(ortho$HumanSymbol), rownames(human_pseudocell)), ] <- human_pseudocell[intersect(as.character(ortho$HumanSymbol), rownames(human_pseudocell)), ]
human_pseudocell_flt <- as.matrix(human_pseudocell_flt)
rownames(human_pseudocell_flt) <- as.character(ortho$FlyGeneID) # change "target" gene id to "ref" gene id

human_label <- DataFrame(row.names = colnames(human_pseudocell_flt), 
                             celltype = unlist(lapply(colnames(human_pseudocell_flt), function (x) unlist(strsplit(as.character(x), split = ' pc'))[1] )) )
human_label$study_id <- 'Human'
head(human_label); nrow(human_label) # 26275


## Ref + Target ###
merged <- cbind(ref_pseudocell_flt, 
                fish_pseudocell_flt, 
                mouse_pseudocell_flt,
                human_pseudocell_flt); dim(merged)
merged_label <- rbind(ref_label, fish_label, mouse_label, human_label)
head(merged_label); nrow(merged_label)
merged_se <- SummarizedExperiment(assays = list(counts=merged), colData = merged_label); merged_se
#saveRDS(merged_se, 'merged_celltype_coregenes.Rds')

