
### Fly-Fish ###
gsva_flytofish <- read.delim('../5.enrichmenttest/1.flytofish/orthologous_allgenes/gsva_flytofish_celltype_v1.txt', row.names = 1, check.names = F)
colnames(gsva_flytofish) <- unlist(lapply(colnames(gsva_flytofish), function (x) paste0(c('Fly|', as.character(x)), collapse = '') ))
rownames(gsva_flytofish) <- unlist(lapply(rownames(gsva_flytofish), function (x) paste0(c('Fish|', as.character(x)), collapse = '') ))
for (ct in colnames(gsva_flytofish)) {
  gsva_flytofish[, ct] <- (gsva_flytofish[, ct] - min(gsva_flytofish[, ct])) / (max(gsva_flytofish[, ct]) - min(gsva_flytofish[, ct]))
}
head(gsva_flytofish)


### Fish-Mouse ###
# MCA
gsva_fishtomouse_mca <- read.delim('../5.enrichmenttest/2.fishtomouse/noNeu_orthologous_allgenes/gsva_fishtomouse_celltype.mca.txt', row.names = 1, check.names = F)
colnames(gsva_fishtomouse_mca) <- unlist(lapply(colnames(gsva_fishtomouse_mca), function (x) paste0(c('Fish|', as.character(x)), collapse = '') ))
rownames(gsva_fishtomouse_mca) <- unlist(lapply(rownames(gsva_fishtomouse_mca), function (x) unlist(strsplit(as.character(x), split = ' \\('))[1] ))
rownames(gsva_fishtomouse_mca) <- unlist(lapply(rownames(gsva_fishtomouse_mca), function (x) paste0(c('Mouse|', as.character(x)), collapse = '') ))
#for (ct in colnames(gsva_fishtomouse_mca)) {
#  gsva_fishtomouse_mca[, ct] <- (gsva_fishtomouse_mca[, ct]-min(gsva_fishtomouse_mca[, ct]))/(max(gsva_fishtomouse_mca[, ct])-min(gsva_fishtomouse_mca[, ct]))
#}
head(gsva_fishtomouse_mca)
# TM 10X
gsva_fishtomouse_tm10x <- read.delim('../5.enrichmenttest/2.fishtomouse/orthologous_allgenes/gsva_fishtomouse_celltype.tm10x.txt', row.names = 1, check.names = F)
colnames(gsva_fishtomouse_tm10x) <- unlist(lapply(colnames(gsva_fishtomouse_tm10x), function (x) paste0(c('Fish|', as.character(x)), collapse = '') ))
rownames(gsva_fishtomouse_tm10x) <- unlist(lapply(rownames(gsva_fishtomouse_tm10x), function (x) unlist(strsplit(as.character(x), split = ' \\('))[1] ))
rownames(gsva_fishtomouse_tm10x) <- unlist(lapply(rownames(gsva_fishtomouse_tm10x), function (x) paste0(c('Mouse|', as.character(x)), collapse = '') ))
#for (ct in colnames(gsva_fishtomouse_tm10x)) {
#  gsva_fishtomouse_tm10x[, ct] <- (gsva_fishtomouse_tm10x[, ct]-min(gsva_fishtomouse_tm10x[, ct]))/(max(gsva_fishtomouse_tm10x[, ct])-min(gsva_fishtomouse_tm10x[, ct]))
#}
head(gsva_fishtomouse_tm10x)
# TM SS2
gsva_fishtomouse_tmss2 <- read.delim('../5.enrichmenttest/2.fishtomouse/orthologous_allgenes/gsva_fishtomouse_celltype.tmss2.txt', row.names = 1, check.names = F)
colnames(gsva_fishtomouse_tmss2) <- unlist(lapply(colnames(gsva_fishtomouse_tmss2), function (x) paste0(c('Fish|', as.character(x)), collapse = '') ))
rownames(gsva_fishtomouse_tmss2) <- unlist(lapply(rownames(gsva_fishtomouse_tmss2), function (x) unlist(strsplit(as.character(x), split = ' \\('))[1] ))
rownames(gsva_fishtomouse_tmss2) <- unlist(lapply(rownames(gsva_fishtomouse_tmss2), function (x) paste0(c('Mouse|', as.character(x)), collapse = '') ))
#for (ct in colnames(gsva_fishtomouse_tmss2)) {
#  gsva_fishtomouse_tmss2[, ct] <- (gsva_fishtomouse_tmss2[, ct]-min(gsva_fishtomouse_tmss2[, ct]))/(max(gsva_fishtomouse_tmss2[, ct])-min(gsva_fishtomouse_tmss2[, ct]))
#}
head(gsva_fishtomouse_tmss2)
# merged 
unique(c(rownames(gsva_fishtomouse_mca), 
         unique(c(rownames(gsva_fishtomouse_tm10x), 
                  rownames(gsva_fishtomouse_tmss2))))) # 13

gsva_fishtomouse <- data.frame(matrix(nrow = 13, ncol = ncol(gsva_fishtomouse_tm10x)))
colnames(gsva_fishtomouse) <- colnames(gsva_fishtomouse_tm10x)
rownames(gsva_fishtomouse) <- c("Mouse|Progenitors", "Mouse|Granulocytopoietic cells", "Mouse|Erythroids", "Mouse|Neutrophils", 
                                "Mouse|Monocytes", "Mouse|Macrophage", "Mouse|DC", "Mouse|pDC", "Mouse|Basophil", 
                                "Mouse|T cells", "Mouse|NK cells", "Mouse|B cells", "Mouse|Plasma cells")

three <- intersect(rownames(gsva_fishtomouse_mca), 
                   intersect(rownames(gsva_fishtomouse_tm10x), rownames(gsva_fishtomouse_tmss2))); three
mca_tm10x <- intersect(rownames(gsva_fishtomouse_mca),
                       setdiff(rownames(gsva_fishtomouse_tm10x), rownames(gsva_fishtomouse_tmss2))); mca_tm10x
mca_tmss2 <- intersect(rownames(gsva_fishtomouse_mca),
                       setdiff(rownames(gsva_fishtomouse_tmss2), rownames(gsva_fishtomouse_tm10x))); mca_tmss2
tm10x_ss2 <- intersect(rownames(gsva_fishtomouse_tm10x),
                       setdiff(rownames(gsva_fishtomouse_tmss2), rownames(gsva_fishtomouse_mca))); tm10x_ss2
mca_unique <- setdiff( rownames(gsva_fishtomouse_mca), c(rownames(gsva_fishtomouse_tm10x), rownames(gsva_fishtomouse_tmss2))); mca_unique

gsva_fishtomouse[three, ] <- (gsva_fishtomouse_mca[three, ] + gsva_fishtomouse_tm10x[three, ] + gsva_fishtomouse_tmss2[three, ])/3
gsva_fishtomouse[mca_tm10x, ] <- (gsva_fishtomouse_mca[mca_tm10x, ] + gsva_fishtomouse_tm10x[mca_tm10x, ])/2
gsva_fishtomouse[mca_tmss2, ] <- (gsva_fishtomouse_mca[mca_tmss2, ] + gsva_fishtomouse_tmss2[mca_tmss2, ])/2
gsva_fishtomouse[tm10x_ss2, ] <- (gsva_fishtomouse_tm10x[tm10x_ss2, ] + gsva_fishtomouse_tmss2[tm10x_ss2, ])/2
gsva_fishtomouse[mca_unique, ] <- gsva_fishtomouse_mca[mca_unique, ]
gsva_fishtomouse

for (ct in colnames(gsva_fishtomouse)) {
  gsva_fishtomouse[, ct] <- (gsva_fishtomouse[, ct]-min(gsva_fishtomouse[, ct]))/(max(gsva_fishtomouse[, ct])-min(gsva_fishtomouse[, ct]))
}
gsva_fishtomouse


### Mouse-Human ###
# HCL
gsva_mousetohuman_hcl <- read.delim('../5.enrichmenttest/3.mousetohuman/noNeu_orthologous_allgenes/gsva_mousetohuman_celltype.hcl.txt', row.names = 1, check.names = F)
colnames(gsva_mousetohuman_hcl) <- unlist(lapply(colnames(gsva_mousetohuman_hcl), function (x) paste0(c('Mouse|', as.character(x)), collapse = '') ))
rownames(gsva_mousetohuman_hcl) <- unlist(lapply(rownames(gsva_mousetohuman_hcl), function (x) unlist(strsplit(as.character(x), split = ' \\('))[1] ))
rownames(gsva_mousetohuman_hcl) <- unlist(lapply(rownames(gsva_mousetohuman_hcl), function (x) paste0(c('Human|', as.character(x)), collapse = '') )); rownames(gsva_mousetohuman_hcl)[11] <- 'Human|Plasma cells'
gsva_mousetohuman_hcl

gsva_mousetohuman_hcl_v2 <- data.frame(matrix(ncol = nrow(gsva_fishtomouse), nrow = nrow(gsva_mousetohuman_hcl)))
colnames(gsva_mousetohuman_hcl_v2) <- rownames(gsva_fishtomouse)
rownames(gsva_mousetohuman_hcl_v2) <- rownames(gsva_mousetohuman_hcl)
gsva_mousetohuman_hcl_v2[, three] <- (gsva_mousetohuman_hcl[, paste(three, ' (MCA)', sep = '')] + 
                                        gsva_mousetohuman_hcl[, paste(three, ' (TM10X)', sep = '')] + gsva_mousetohuman_hcl[, paste(three, ' (TMSS2)', sep = '')])/3
gsva_mousetohuman_hcl_v2[, mca_tm10x] <- (gsva_mousetohuman_hcl[, paste(mca_tm10x, ' (MCA)', sep = '')] + gsva_mousetohuman_hcl[, paste(mca_tm10x, ' (TM10X)', sep = '')])/2
gsva_mousetohuman_hcl_v2[, mca_tmss2] <- (gsva_mousetohuman_hcl[, paste(mca_tmss2, ' (MCA)', sep = '')] + gsva_mousetohuman_hcl[, paste(mca_tmss2, ' (TMSS2)', sep = '')])/2
gsva_mousetohuman_hcl_v2[, tm10x_ss2] <- (gsva_mousetohuman_hcl[, paste(tm10x_ss2, ' (TM10X)', sep = '')] + gsva_mousetohuman_hcl[, paste(tm10x_ss2, ' (TMSS2)', sep = '')])/2
gsva_mousetohuman_hcl_v2[, mca_unique] <- gsva_mousetohuman_hcl[, paste(mca_unique, ' (MCA)', sep = '')]

#for (ct in colnames(gsva_mousetohuman_hcl_v2)) {
#  gsva_mousetohuman_hcl_v2[, ct] <- (gsva_mousetohuman_hcl_v2[, ct]-min(gsva_mousetohuman_hcl_v2[, ct]))/(max(gsva_mousetohuman_hcl_v2[, ct])-min(gsva_mousetohuman_hcl_v2[, ct]))
#}
head(gsva_mousetohuman_hcl_v2)

# HCA
gsva_mousetohuman_hca <- read.delim('../5.enrichmenttest/3.mousetohuman/noNeu_orthologous_allgenes/gsva_mousetohuman_celltype.hca.txt', row.names = 1, check.names = F)
colnames(gsva_mousetohuman_hca) <- unlist(lapply(colnames(gsva_mousetohuman_hca), function (x) paste0(c('Mouse|', as.character(x)), collapse = '') ))
rownames(gsva_mousetohuman_hca) <- unlist(lapply(rownames(gsva_mousetohuman_hca), function (x) unlist(strsplit(as.character(x), split = ' \\('))[1] ))
rownames(gsva_mousetohuman_hca) <- unlist(lapply(rownames(gsva_mousetohuman_hca), function (x) paste0(c('Human|', as.character(x)), collapse = '') ))
gsva_mousetohuman_hca

gsva_mousetohuman_hca_v2 <- data.frame(matrix(ncol = nrow(gsva_fishtomouse), nrow = nrow(gsva_mousetohuman_hca)))
colnames(gsva_mousetohuman_hca_v2) <- rownames(gsva_fishtomouse)
rownames(gsva_mousetohuman_hca_v2) <- rownames(gsva_mousetohuman_hca)
gsva_mousetohuman_hca_v2[, three] <- (gsva_mousetohuman_hca[, paste(three, ' (MCA)', sep = '')] + 
                                        gsva_mousetohuman_hca[, paste(three, ' (TM10X)', sep = '')] + gsva_mousetohuman_hca[, paste(three, ' (TMSS2)', sep = '')])/3
gsva_mousetohuman_hca_v2[, mca_tm10x] <- (gsva_mousetohuman_hca[, paste(mca_tm10x, ' (MCA)', sep = '')] + gsva_mousetohuman_hca[, paste(mca_tm10x, ' (TM10X)', sep = '')])/2
gsva_mousetohuman_hca_v2[, mca_tmss2] <- (gsva_mousetohuman_hca[, paste(mca_tmss2, ' (MCA)', sep = '')] + gsva_mousetohuman_hca[, paste(mca_tmss2, ' (TMSS2)', sep = '')])/2
gsva_mousetohuman_hca_v2[, tm10x_ss2] <- (gsva_mousetohuman_hca[, paste(tm10x_ss2, ' (TM10X)', sep = '')] + gsva_mousetohuman_hca[, paste(tm10x_ss2, ' (TMSS2)', sep = '')])/2
gsva_mousetohuman_hca_v2[, mca_unique] <- gsva_mousetohuman_hca[, paste(mca_unique, ' (MCA)', sep = '')]
gsva_mousetohuman_hca_v2

#for (ct in colnames(gsva_mousetohuman_hca_v2)) {
#  gsva_mousetohuman_hca_v2[, ct] <- (gsva_mousetohuman_hca_v2[, ct]-min(gsva_mousetohuman_hca_v2[, ct]))/(max(gsva_mousetohuman_hca_v2[, ct])-min(gsva_mousetohuman_hca_v2[, ct]))
#}
head(gsva_mousetohuman_hca_v2)

# merged 
unique(c( rownames(gsva_mousetohuman_hcl_v2), rownames(gsva_mousetohuman_hca_v2) )) # 17

gsva_mousetohuman <- data.frame(matrix(nrow = 17, ncol = ncol(gsva_mousetohuman_hcl_v2)))
colnames(gsva_mousetohuman) <- colnames(gsva_mousetohuman_hcl_v2)
rownames(gsva_mousetohuman) <- c("Human|Progenitors", "Human|Granulocyte progenitor", "Human|Erythroids", "Human|Neutrophils", "Human|Monocytes", "Human|Macrophage", "Human|DC", "Human|pDC", 
                                 "Human|T cells", "Human|CD8 T cells", "Human|NK cells", "Human|pre-PC", "Human|pro-B cells", "Human|pre-B cells", "Human|B cells", "Human|Plasma cells", "Human|Platelet")

both <- intersect(rownames(gsva_mousetohuman_hcl_v2), rownames(gsva_mousetohuman_hca_v2)); both
hcl_unique <- setdiff(rownames(gsva_mousetohuman_hcl_v2), rownames(gsva_mousetohuman_hca_v2)); hcl_unique
hca_unique <- setdiff(rownames(gsva_mousetohuman_hca_v2), rownames(gsva_mousetohuman_hcl_v2)); hca_unique

gsva_mousetohuman[both, ] <- (gsva_mousetohuman_hcl_v2[both, ] + gsva_mousetohuman_hca_v2[both, ])/2
gsva_mousetohuman[hcl_unique, ] <- gsva_mousetohuman_hcl_v2[hcl_unique, ]
gsva_mousetohuman[hca_unique, ] <- gsva_mousetohuman_hca_v2[hca_unique, ]
gsva_mousetohuman

for (ct in colnames(gsva_mousetohuman)) {
  gsva_mousetohuman[, ct] <- (gsva_mousetohuman[, ct]-min(gsva_mousetohuman[, ct]))/(max(gsva_mousetohuman[, ct])-min(gsva_mousetohuman[, ct]))
}
gsva_mousetohuman


meta_flytofish <- read.delim('../6.MetaNeighbor/1.flytofish/celltype/MetaNeighbor_flytofish_orthologous.celltype_top_hits.txt')
for (i in c(1:nrow(meta_flytofish)) ) {
  ct1 <- as.character(meta_flytofish[i, 1])
  ct2 <- as.character(meta_flytofish[i, 2])
  gsvascore <- c(gsva_flytofish[ct1, ct2], gsva_flytofish[ct2, ct1])
  meta_flytofish[i, 'avgscore'] <- mean(c(meta_flytofish[i,3], gsvascore))
}
head(meta_flytofish)

meta_fishtomouse <- read.delim('../6.MetaNeighbor/2.fishtomouse_noNeu/celltype_merged/MetaNeighbor_flytofish_orthologous.celltype_top_hits.txt')
for (i in c(1:nrow(meta_fishtomouse)) ) {
  ct1 <- as.character(meta_fishtomouse[i, 1])
  ct2 <- as.character(meta_fishtomouse[i, 2])
  gsvascore <- c(gsva_fishtomouse[ct1, ct2], gsva_fishtomouse[ct2, ct1])
  meta_fishtomouse[i, 'avgscore'] <- mean(c(meta_fishtomouse[i,3], gsvascore))
}
head(meta_fishtomouse)

meta_mousetohuman <- read.delim('../6.MetaNeighbor/3.mousetohuman_noNeu/celltype_merged/MetaNeighbor_orthologous.celltype_top_hits.txt')
for (i in c(1:nrow(meta_mousetohuman)) ) {
  ct1 <- as.character(meta_mousetohuman[i, 1])
  ct2 <- as.character(meta_mousetohuman[i, 2])
  gsvascore <- c(gsva_mousetohuman[ct1, ct2], gsva_mousetohuman[ct2, ct1])
  meta_mousetohuman[i, 'avgscore'] <- mean(c(meta_mousetohuman[i,3], gsvascore))
}
head(meta_mousetohuman)

###
merged <- rbind(meta_flytofish, meta_fishtomouse, meta_mousetohuman)
colnames(merged) <- c('Celltype_1', 'Celltype_2', 'Mean_AUROC', 'Match_type', 'Avgscore')
merged$Celltype_1 <- as.character(merged$Celltype_1)
merged$Celltype_2 <- as.character(merged$Celltype_2)
merged$Match_type <- as.character(merged$Match_type)

merged <- subset(merged, Avgscore >= 0.80 | Match_type == 'Reciprocal_top_hit')
rownames(merged) <- c(1:nrow(merged))
merged

merged[2, c(1,2)] <- merged[2, c(2,1)]
merged[c(12,13), c(1,2)] <- merged[c(12,13), c(2,1)]
merged[15, c(1,2)] <- merged[15, c(2,1)]
merged[21, c(1,2)] <- merged[21, c(2,1)]
merged[24, c(1,2)] <- merged[24, c(2,1)]
merged[c(27,28,30:33), c(1,2)] <- merged[c(27,28,30:33), c(2,1)]
merged

merged$hierarchy <- c('fzf','fzf','fzf','fzf','fzf', 'fzf', 
                      'fm', 'fm', 'fm', 'fm', 'fm', 'fm', 'fm', 'fm', 'fm', 
                      'mh','mh','mh','mh','mh','mh','mh', 'mh','mh','mh','mh','mh','mh','mh', 'mh','mh','mh','mh')
merged$Celltype_1_group <- c('Fly', 'Fly', 'Fly', 'Fly', 'Fly', 'Fly', 
                             'Fish', 'Fish', 'Fish', 'Fish', 'Fish', 'Fish', 'Fish', 'Fish', 'Fish', 
                             'Mouse', 'Mouse', 'Mouse', 'Mouse', 'Mouse', 'Mouse', 'Mouse', 'Mouse', 'Mouse', 'Mouse', 'Mouse', 'Mouse', 'Mouse', 'Mouse', 'Mouse', 'Mouse', 'Mouse', 'Mouse')
merged$Celltype_2_group <- c('Fish', 'Fish', 'Fish', 'Fish', 'Fish', 'Fish', 
                             'Mouse', 'Mouse', 'Mouse', 'Mouse', 'Mouse', 'Mouse', 'Mouse', 'Mouse', 'Mouse', 
                             'Human', 'Human', 'Human', 'Human', 'Human', 'Human', 'Human', 'Human', 'Human', 'Human', 'Human', 'Human', 'Human', 'Human', 'Human', 'Human', 'Human', 'Human')
merged
#write.table(merged, 'heirarchy_noNeu.txt', quote = F, sep = '\t', row.names = F, col.names = T)


evol <- data.frame(matrix(ncol = ncol(merged), nrow = ncol(gsva_flytofish)*nrow(gsva_flytofish) + ncol(gsva_fishtomouse)*nrow(gsva_fishtomouse) + nrow(gsva_fishtomouse)*nrow(gsva_mousetohuman) ))
dim(evol)
colnames(evol) <- colnames(merged)

evol$Celltype_1 <- c( rep(colnames(gsva_flytofish), each = nrow(gsva_flytofish)), # 48
                      rep(colnames(gsva_fishtomouse), each = nrow(gsva_fishtomouse)), # 78
                      rep(rownames(gsva_fishtomouse), each = nrow(gsva_mousetohuman))) # 221

evol$Celltype_2 <- c( rep(rownames(gsva_flytofish), times = ncol(gsva_flytofish)), # 48
                      rep(rownames(gsva_fishtomouse), times = ncol(gsva_fishtomouse)), # 78
                      rep(rownames(gsva_mousetohuman), times = ncol(gsva_mousetohuman))) # 221
evol$Mean_AUROC <- 'NS'
evol$Match_type <- 'NS'
evol$Avgscore <- 0
evol$hierarchy <- c( rep('fzf', each = ncol(gsva_flytofish)*nrow(gsva_flytofish)),
                     rep('fm', each = ncol(gsva_fishtomouse)*nrow(gsva_fishtomouse)),
                     rep('mh', each = nrow(gsva_fishtomouse)*nrow(gsva_mousetohuman)))
evol$Celltype_1_group <- c( rep('Fly', each = ncol(gsva_flytofish)*nrow(gsva_flytofish)),
                            rep('Fish', each = ncol(gsva_fishtomouse)*nrow(gsva_fishtomouse)),
                            rep('Mouse', each = nrow(gsva_fishtomouse)*nrow(gsva_mousetohuman)))
evol$Celltype_2_group <- c( rep('Fish', each = ncol(gsva_flytofish)*nrow(gsva_flytofish)),
                            rep('Mouse', each = ncol(gsva_fishtomouse)*nrow(gsva_fishtomouse)),
                            rep('Human', each = nrow(gsva_fishtomouse)*nrow(gsva_mousetohuman)))

for (conn in rownames(merged)) {
  evol[evol$Celltype_1 == merged[conn, 'Celltype_1'] & evol$Celltype_2 == merged[conn, 'Celltype_2'], ] <- merged[conn, ]
}
head(evol)
#write.table(evol, 'heirarchy_noNeu_cytoscape.txt', quote = F, sep = '\t', row.names = F, col.names = T)



merged_flt <- data.frame(matrix(nrow = 0, ncol = ncol(merged)))
colnames(merged_flt) <- colnames(merged)
for (unique_type in unique(merged$Celltype_1)) {
  tmp <- subset(merged, Celltype_1 == unique_type)
  if (nrow(tmp) == 1) {
    merged_flt <- rbind(merged_flt, tmp)
  } else {
    merged_flt <- rbind(merged_flt, subset(tmp, Avgscore == max(tmp$Avgscore)) )
  }
}
rownames(merged_flt) <- c(1:nrow(merged_flt))
merged_flt
#write.table(merged_flt, 'heirarchy_noNeu.flt.txt', quote = F, sep = '\t', row.names = F, col.names = T)


evol <- data.frame(matrix(ncol = ncol(merged_flt), nrow = ncol(gsva_flytofish)*nrow(gsva_flytofish) + ncol(gsva_fishtomouse)*nrow(gsva_fishtomouse) + nrow(gsva_fishtomouse)*nrow(gsva_mousetohuman) ))
dim(evol)
colnames(evol) <- colnames(merged_flt)

evol$Celltype_1 <- c( rep(colnames(gsva_flytofish), each = nrow(gsva_flytofish)), # 48
                      rep(colnames(gsva_fishtomouse), each = nrow(gsva_fishtomouse)), # 78
                      rep(rownames(gsva_fishtomouse), each = nrow(gsva_mousetohuman))) # 221

evol$Celltype_2 <- c( rep(rownames(gsva_flytofish), times = ncol(gsva_flytofish)), # 48
                      rep(rownames(gsva_fishtomouse), times = ncol(gsva_fishtomouse)), # 78
                      rep(rownames(gsva_mousetohuman), times = ncol(gsva_mousetohuman))) # 221
evol$Mean_AUROC <- 'NS'
evol$Match_type <- 'NS'
evol$Avgscore <- 0
evol$hierarchy <- c( rep('fzf', each = ncol(gsva_flytofish)*nrow(gsva_flytofish)),
                     rep('fm', each = ncol(gsva_fishtomouse)*nrow(gsva_fishtomouse)),
                     rep('mh', each = nrow(gsva_fishtomouse)*nrow(gsva_mousetohuman)))
evol$Celltype_1_group <- c( rep('Fly', each = ncol(gsva_flytofish)*nrow(gsva_flytofish)),
                            rep('Fish', each = ncol(gsva_fishtomouse)*nrow(gsva_fishtomouse)),
                            rep('Mouse', each = nrow(gsva_fishtomouse)*nrow(gsva_mousetohuman)))
evol$Celltype_2_group <- c( rep('Fish', each = ncol(gsva_flytofish)*nrow(gsva_flytofish)),
                            rep('Mouse', each = ncol(gsva_fishtomouse)*nrow(gsva_fishtomouse)),
                            rep('Human', each = nrow(gsva_fishtomouse)*nrow(gsva_mousetohuman)))

for (conn in rownames(merged_flt)) {
  evol[evol$Celltype_1 == merged_flt[conn, 'Celltype_1'] & evol$Celltype_2 == merged_flt[conn, 'Celltype_2'], ] <- merged_flt[conn, ]
}
head(evol)
#write.table(evol, 'heirarchy_noNeu_cytoscape.flt.txt', quote = F, sep = '\t', row.names = F, col.names = T)


