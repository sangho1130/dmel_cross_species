
### Fly-Fish ###
gsva_flytofish <- read.delim('../5.enrichmenttest/1.flytofish/orthologous_allgenes/gsva_flytofish_celltype_v1.txt', row.names = 1, check.names = F)
for (ct in colnames(gsva_flytofish)) {
  gsva_flytofish[, ct] <- (gsva_flytofish[, ct] - min(gsva_flytofish[, ct])) / (max(gsva_flytofish[, ct]) - min(gsva_flytofish[, ct]))
}
head(gsva_flytofish)


### Fish-Mouse ###
gsva_fishtomouse <- read.delim('../5.enrichmenttest/2.fishtomouse_v2/orthologous_allgenes/gsva_fishtomouse_celltype.txt', row.names = 1, check.names = F)
for (ct in colnames(gsva_fishtomouse)) {
  gsva_fishtomouse[, ct] <- (gsva_fishtomouse[, ct]-min(gsva_fishtomouse[, ct]))/(max(gsva_fishtomouse[, ct])-min(gsva_fishtomouse[, ct]))
}
head(gsva_fishtomouse)

### Mouse-Human ###
gsva_mousetohuman <- read.delim('../5.enrichmenttest/3.mousetohuman/noPlt_orthologous_allgenes/gsva_mousetohuman_celltype.txt', row.names = 1, check.names = F)
for (ct in colnames(gsva_mousetohuman)) {
  gsva_mousetohuman[, ct] <- (gsva_mousetohuman[, ct]-min(gsva_mousetohuman[, ct]))/(max(gsva_mousetohuman[, ct])-min(gsva_mousetohuman[, ct]))
}
head(gsva_mousetohuman)


### Fly-Fish ###
meta_flytofish <- read.delim('../6.MetaNeighbor/1.flytofish/celltype/MetaNeighbor_flytofish_orthologous.celltype_top_hits.txt')
for (i in c(1:nrow(meta_flytofish)) ) {
  ct1 <- unlist(strsplit(as.character(meta_flytofish[i, 1]), split = '\\|'))[2]
  ct2 <- unlist(strsplit(as.character(meta_flytofish[i, 2]), split = '\\|'))[2]
  
  if (ct1 %in% rownames(gsva_flytofish)) {
    gsvascore <- gsva_flytofish[ct1, ct2]
  } else {
    gsvascore <- gsva_flytofish[ct2, ct1]
  }
  
  meta_flytofish[i, 'avgscore'] <- mean(c(meta_flytofish[i,3], gsvascore))
}
head(meta_flytofish)

### Fish-Mouse ###
meta_fishtomouse <- read.delim('../6.MetaNeighbor/2.fishtomouse_v2/celltype/MetaNeighbor_fishtomouse_orthologous.celltype_top_hits.txt')
for (i in c(1:nrow(meta_fishtomouse)) ) {
  ct1 <- unlist(strsplit(as.character(meta_fishtomouse[i, 1]), split = '\\|'))[2]
  ct2 <- unlist(strsplit(as.character(meta_fishtomouse[i, 2]), split = '\\|'))[2]
  
  if (ct1 %in% rownames(gsva_fishtomouse)) {
    gsvascore <- gsva_fishtomouse[ct1, ct2]
  } else {
    gsvascore <- gsva_fishtomouse[ct2, ct1]
  }
  
  meta_fishtomouse[i, 'avgscore'] <- mean(c(meta_fishtomouse[i,3], gsvascore))
}
head(meta_fishtomouse)

### Mouse-Human ###
meta_mousetohuman <- read.delim('../6.MetaNeighbor/3.mousetohuman_v2/celltype_noPlt/MetaNeighbor_orthologous.celltype_top_hits.txt')
for (i in c(1:nrow(meta_mousetohuman)) ) {
  ct1 <- unlist(strsplit(as.character(meta_mousetohuman[i, 1]), split = '\\|'))[2]
  ct2 <- unlist(strsplit(as.character(meta_mousetohuman[i, 2]), split = '\\|'))[2]
  
  if (ct1 %in% rownames(gsva_mousetohuman)) {
    gsvascore <- gsva_mousetohuman[ct1, ct2]
  } else {
    gsvascore <- gsva_mousetohuman[ct2, ct1]
  }
  
  meta_mousetohuman[i, 'avgscore'] <- mean(c(meta_mousetohuman[i,3], gsvascore))
}
head(meta_mousetohuman)

###
merged <- rbind(meta_flytofish, meta_fishtomouse, meta_mousetohuman)
colnames(merged) <- c('Celltype_1', 'Celltype_2', 'Mean_AUROC', 'Match_type', 'Avgscore')
merged$Celltype_1 <- as.character(merged$Celltype_1)
merged$Celltype_2 <- as.character(merged$Celltype_2)
merged$Match_type <- as.character(merged$Match_type)
merged

merged <- subset(merged, Avgscore >= 0.80 | Match_type == 'Reciprocal_top_hit')
rownames(merged) <- c(1:nrow(merged))
merged # 31

merged[c(2), c(1, 2)] <- merged[c(2), c(2, 1)]
merged[c(12, 13), c(1, 2)] <- merged[c(12, 13), c(2, 1)]
merged[c(18, 19, 21, 24, 26, 27, 30), c(1, 2)] <- merged[c(18, 19, 21, 24, 26, 27, 30), c(2, 1)]
merged

merged$hierarchy <- unlist(lapply(c(1:nrow(merged)), function (x) paste( unlist(strsplit(merged$Celltype_1[x], split = '|'))[1], unlist(strsplit(merged$Celltype_2[x], split = '|'))[1], sep = '') ))
merged$Celltype_1_group <- unlist(lapply(c(1:nrow(merged)), function (x) unlist(strsplit(merged$Celltype_1[x], split = '\\|'))[1]))
merged$Celltype_2_group <- unlist(lapply(c(1:nrow(merged)), function (x) unlist(strsplit(merged$Celltype_2[x], split = '\\|'))[1]))
merged$Celltype_1 <- unlist(lapply(c(1:nrow(merged)), function (x) unlist(strsplit(merged$Celltype_1[x], split = '\\|'))[2]))
merged$Celltype_2 <- unlist(lapply(c(1:nrow(merged)), function (x) unlist(strsplit(merged$Celltype_2[x], split = '\\|'))[2]))
merged
#write.table(merged, 'heirarchy_merged.txt', quote = F, sep = '\t', row.names = F, col.names = T)


evol <- data.frame(matrix(ncol = ncol(merged), nrow = ncol(gsva_flytofish)*nrow(gsva_flytofish) + ncol(gsva_fishtomouse)*nrow(gsva_fishtomouse) + nrow(gsva_fishtomouse)*nrow(gsva_mousetohuman) ))
dim(evol)
colnames(evol) <- colnames(merged)


evol$Celltype_1 <- c( rep(colnames(gsva_flytofish), each = nrow(gsva_flytofish)), # 48
                      rep(colnames(gsva_fishtomouse), each = nrow(gsva_fishtomouse)), # 78
                      rep(rownames(gsva_fishtomouse), each = nrow(gsva_mousetohuman))) # 208

evol$Celltype_2 <- c( rep(rownames(gsva_flytofish), times = ncol(gsva_flytofish)), # 48
                      rep(rownames(gsva_fishtomouse), times = ncol(gsva_fishtomouse)), # 78
                      rep(rownames(gsva_mousetohuman), times = ncol(gsva_mousetohuman))) # 208
evol$Mean_AUROC <- 'NS'
evol$Match_type <- 'NS'
evol$Avgscore <- 0
evol$hierarchy <- c( rep('FF', each = ncol(gsva_flytofish)*nrow(gsva_flytofish)),
                     rep('FM', each = ncol(gsva_fishtomouse)*nrow(gsva_fishtomouse)),
                     rep('MH', each = nrow(gsva_fishtomouse)*nrow(gsva_mousetohuman)))
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
#write.table(evol, 'heirarchy_merged_cytoscape.txt', quote = F, sep = '\t', row.names = F, col.names = T)



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
#write.table(merged_flt, 'heirarchy_merged_cytoscape.flt.txt', quote = F, sep = '\t', row.names = F, col.names = T)


evol <- data.frame(matrix(ncol = ncol(merged_flt), nrow = ncol(gsva_flytofish)*nrow(gsva_flytofish) + ncol(gsva_fishtomouse)*nrow(gsva_fishtomouse) + nrow(gsva_fishtomouse)*nrow(gsva_mousetohuman) ))
dim(evol)
colnames(evol) <- colnames(merged_flt)

evol$Celltype_1 <- c( rep(colnames(gsva_flytofish), each = nrow(gsva_flytofish)), # 48
                      rep(colnames(gsva_fishtomouse), each = nrow(gsva_fishtomouse)), # 78
                      rep(rownames(gsva_fishtomouse), each = nrow(gsva_mousetohuman))) # 208

evol$Celltype_2 <- c( rep(rownames(gsva_flytofish), times = ncol(gsva_flytofish)), # 48
                      rep(rownames(gsva_fishtomouse), times = ncol(gsva_fishtomouse)), # 78
                      rep(rownames(gsva_mousetohuman), times = ncol(gsva_mousetohuman))) # 208
evol$Mean_AUROC <- 'NS'
evol$Match_type <- 'NS'
evol$Avgscore <- 0
evol$hierarchy <- c( rep('FF', each = ncol(gsva_flytofish)*nrow(gsva_flytofish)),
                     rep('FM', each = ncol(gsva_fishtomouse)*nrow(gsva_fishtomouse)),
                     rep('MH', each = nrow(gsva_fishtomouse)*nrow(gsva_mousetohuman)))
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
#write.table(evol, 'heirarchy_merged_cytoscape.flt.txt', quote = F, sep = '\t', row.names = F, col.names = T)


