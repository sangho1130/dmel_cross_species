
### Fly-Fish ###
gsva_flytofish <- read.delim('../5.enrichmenttest/1.flytofish/orthologous_allgenes/gsva_flytofish_celltype_v1.txt', row.names = 1, check.names = F)
colnames(gsva_flytofish) <- unlist(lapply(colnames(gsva_flytofish), function (x) paste0(c('Fly|', as.character(x)), collapse = '') ))
rownames(gsva_flytofish) <- unlist(lapply(rownames(gsva_flytofish), function (x) paste0(c('Fish|', as.character(x)), collapse = '') ))
for (ct in colnames(gsva_flytofish)) {
  gsva_flytofish[, ct] <- (gsva_flytofish[, ct] - min(gsva_flytofish[, ct])) / (max(gsva_flytofish[, ct]) - min(gsva_flytofish[, ct]))
}
gsva_flytofish
#################


### Fly-Mouse ###
gsva_flytomouse <- read.delim('../5.enrichmenttest/1.flytomouse_v2/orthologous_allgenes/gsva_flytomouse_celltype.txt', row.names = 1, check.names = F)
colnames(gsva_flytomouse) <- unlist(lapply(colnames(gsva_flytomouse), function (x) paste0(c('Fly|', as.character(x)), collapse = '') ))
rownames(gsva_flytomouse) <- unlist(lapply(rownames(gsva_flytomouse), function (x) unlist(strsplit(as.character(x), split = ' \\('))[1] ))
rownames(gsva_flytomouse) <- unlist(lapply(rownames(gsva_flytomouse), function (x) paste0(c('Mouse|', as.character(x)), collapse = '') ))
head(gsva_flytomouse)

for (ct in colnames(gsva_flytomouse)) {
  gsva_flytomouse[, ct] <- (gsva_flytomouse[, ct]-min(gsva_flytomouse[, ct]))/(max(gsva_flytomouse[, ct])-min(gsva_flytomouse[, ct]))
}
gsva_flytomouse
#################


### Fly-Human ###
gsva_flytohuman <- read.delim('../5.enrichmenttest/1.flytohuman/noPlt_orthologous_allgenes/gsva_flytohuman_celltype.txt', row.names = 1, check.names = F)
colnames(gsva_flytohuman) <- unlist(lapply(colnames(gsva_flytohuman), function (x) paste0(c('Fly|', as.character(x)), collapse = '') ))
rownames(gsva_flytohuman) <- unlist(lapply(rownames(gsva_flytohuman), function (x) unlist(strsplit(as.character(x), split = ' \\('))[1] ))
rownames(gsva_flytohuman) <- unlist(lapply(rownames(gsva_flytohuman), function (x) paste0(c('Human|', as.character(x)), collapse = '') ))
head(gsva_flytohuman)

for (ct in colnames(gsva_flytohuman)) {
  gsva_flytohuman[, ct] <- (gsva_flytohuman[, ct]-min(gsva_flytohuman[, ct]))/(max(gsva_flytohuman[, ct])-min(gsva_flytohuman[, ct]))
}
gsva_flytohuman
##################


### Fish-Mouse ###
gsva_fishtomouse <- read.delim('../5.enrichmenttest/2.fishtomouse_v2/orthologous_allgenes/gsva_fishtomouse_celltype.txt', row.names = 1, check.names = F)
colnames(gsva_fishtomouse) <- unlist(lapply(colnames(gsva_fishtomouse), function (x) paste0(c('Fish|', as.character(x)), collapse = '') ))
rownames(gsva_fishtomouse) <- unlist(lapply(rownames(gsva_fishtomouse), function (x) unlist(strsplit(as.character(x), split = ' \\('))[1] ))
rownames(gsva_fishtomouse) <- unlist(lapply(rownames(gsva_fishtomouse), function (x) paste0(c('Mouse|', as.character(x)), collapse = '') ))
head(gsva_fishtomouse)

for (ct in colnames(gsva_fishtomouse)) {
  gsva_fishtomouse[, ct] <- (gsva_fishtomouse[, ct]-min(gsva_fishtomouse[, ct]))/(max(gsva_fishtomouse[, ct])-min(gsva_fishtomouse[, ct]))
}
gsva_fishtomouse
##################


### Fish-Human ###
gsva_fishtohuman <- read.delim('../5.enrichmenttest/2.fishtohuman/noPlt_orthologous_allgenes/gsva_fishtohuman_celltype.txt', row.names = 1, check.names = F)
colnames(gsva_fishtohuman) <- unlist(lapply(colnames(gsva_fishtohuman), function (x) paste0(c('Fish|', as.character(x)), collapse = '') ))
rownames(gsva_fishtohuman) <- unlist(lapply(rownames(gsva_fishtohuman), function (x) unlist(strsplit(as.character(x), split = ' \\('))[1] ))
rownames(gsva_fishtohuman) <- unlist(lapply(rownames(gsva_fishtohuman), function (x) paste0(c('Human|', as.character(x)), collapse = '') ))
head(gsva_fishtohuman)

for (ct in colnames(gsva_fishtohuman)) {
  gsva_fishtohuman[, ct] <- (gsva_fishtohuman[, ct]-min(gsva_fishtohuman[, ct]))/(max(gsva_fishtohuman[, ct])-min(gsva_fishtohuman[, ct]))
}
gsva_fishtohuman
###################


### Mouse-Human ###
gsva_mousetohuman <- read.delim('../5.enrichmenttest/3.mousetohuman/noPlt_orthologous_allgenes/gsva_mousetohuman_celltype.txt', row.names = 1, check.names = F)
colnames(gsva_mousetohuman) <- unlist(lapply(colnames(gsva_mousetohuman), function (x) paste0(c('Mouse|', as.character(x)), collapse = '') ))
rownames(gsva_mousetohuman) <- unlist(lapply(rownames(gsva_mousetohuman), function (x) unlist(strsplit(as.character(x), split = ' \\('))[1] ))
rownames(gsva_mousetohuman) <- unlist(lapply(rownames(gsva_mousetohuman), function (x) paste0(c('Human|', as.character(x)), collapse = '') ))
gsva_mousetohuman

for (ct in colnames(gsva_mousetohuman)) {
  gsva_mousetohuman[, ct] <- (gsva_mousetohuman[, ct]-min(gsva_mousetohuman[, ct]))/(max(gsva_mousetohuman[, ct])-min(gsva_mousetohuman[, ct]))
}
gsva_mousetohuman
####################


### MetaNeighbor ###
### Fly ###
meta_flytofish <- read.delim('../6.MetaNeighbor/1.flytofish/celltype/MetaNeighbor_flytofish_orthologous.celltype_top_hits.txt')
for (i in c(1:nrow(meta_flytofish)) ) {
  ct1 <- as.character(meta_flytofish[i, 1])
  ct2 <- as.character(meta_flytofish[i, 2])
  gsvascore <- c(gsva_flytofish[ct1, ct2], gsva_flytofish[ct2, ct1])
  meta_flytofish[i, 'avgscore'] <- mean(c(meta_flytofish[i,3], gsvascore))
}
head(meta_flytofish)


meta_flytomouse <- read.delim('../6.MetaNeighbor/1.flytomouse_v2/celltype/MetaNeighbor_flytomouse_orthologous.celltype_top_hits.txt')
for (i in c(1:nrow(meta_flytomouse)) ) {
  ct1 <- as.character(meta_flytomouse[i, 1])
  ct2 <- as.character(meta_flytomouse[i, 2])
  gsvascore <- c(gsva_flytomouse[ct1, ct2], gsva_flytomouse[ct2, ct1])
  meta_flytomouse[i, 'avgscore'] <- mean(c(meta_flytomouse[i,3], gsvascore))
}
head(meta_flytomouse)


meta_flytohuman <- read.delim('../6.MetaNeighbor/1.flytohuman/celltype_noPlt/MetaNeighbor_flytohuman_orthologous.celltype_top_hits.txt')
for (i in c(1:nrow(meta_flytohuman)) ) {
  ct1 <- as.character(meta_flytohuman[i, 1])
  ct2 <- as.character(meta_flytohuman[i, 2])
  gsvascore <- c(gsva_flytohuman[ct1, ct2], gsva_flytohuman[ct2, ct1])
  meta_flytohuman[i, 'avgscore'] <- mean(c(meta_flytohuman[i,3], gsvascore))
}
head(meta_flytohuman)


### Fish ###
meta_fishtomouse <- read.delim('../6.MetaNeighbor/2.fishtomouse_v2/celltype/MetaNeighbor_fishtomouse_orthologous.celltype_top_hits.txt')
for (i in c(1:nrow(meta_fishtomouse)) ) {
  ct1 <- as.character(meta_fishtomouse[i, 1])
  ct2 <- as.character(meta_fishtomouse[i, 2])
  gsvascore <- c(gsva_fishtomouse[ct1, ct2], gsva_fishtomouse[ct2, ct1])
  meta_fishtomouse[i, 'avgscore'] <- mean(c(meta_fishtomouse[i,3], gsvascore))
}
head(meta_fishtomouse)

meta_fishtohuman <- read.delim('../6.MetaNeighbor/2.fishtohuman/celltype_noPlt/MetaNeighbor_fishtohuman_orthologous.celltype_top_hits.txt')
for (i in c(1:nrow(meta_fishtohuman)) ) {
  ct1 <- as.character(meta_fishtohuman[i, 1])
  ct2 <- as.character(meta_fishtohuman[i, 2])
  gsvascore <- c(gsva_fishtohuman[ct1, ct2], gsva_fishtohuman[ct2, ct1])
  meta_fishtohuman[i, 'avgscore'] <- mean(c(meta_fishtohuman[i,3], gsvascore))
}
head(meta_fishtohuman)


### Mouse ###
meta_mousetohuman <- read.delim('../6.MetaNeighbor/3.mousetohuman_v2/celltype_noPlt/MetaNeighbor_orthologous.celltype_top_hits.txt')
for (i in c(1:nrow(meta_mousetohuman)) ) {
  ct1 <- as.character(meta_mousetohuman[i, 1])
  ct2 <- as.character(meta_mousetohuman[i, 2])
  gsvascore <- c(gsva_mousetohuman[ct1, ct2], gsva_mousetohuman[ct2, ct1])
  meta_mousetohuman[i, 'avgscore'] <- mean(c(meta_mousetohuman[i,3], gsvascore))
}
head(meta_mousetohuman)


###
merged <- rbind(meta_flytofish, meta_flytomouse, meta_flytohuman, meta_fishtomouse, meta_fishtohuman, meta_mousetohuman)
colnames(merged) <- c('Celltype_1', 'Celltype_2', 'Mean_AUROC', 'Match_type', 'Avgscore')
merged$Celltype_1 <- as.character(merged$Celltype_1)
merged$Celltype_2 <- as.character(merged$Celltype_2)
merged$Match_type <- as.character(merged$Match_type)

merged <- subset(merged, Avgscore >= 0.80 | Match_type == 'Reciprocal_top_hit')
rownames(merged) <- c(1:nrow(merged))
merged

merged[2, c(1,2)] <- merged[2, c(2,1)]
merged[c(13), c(1,2)] <- merged[c(13), c(2,1)]
merged[c(19, 21, 22), c(1,2)] <- merged[c(19, 21, 22), c(2,1)]
merged[c(28, 29), c(1,2)] <- merged[c(28, 29), c(2,1)]
merged[c(32, 34:37, 41), c(1,2)] <- merged[c(32, 34:37, 41), c(2,1)]
merged[c(46, 47, 49, 52, 55, 56, 59), c(1,2)] <- merged[c(46, 47, 49, 52, 55, 56, 59), c(2,1)]
merged

merged$hierarchy <- paste(unlist(lapply(merged$Celltype_1, function (x) unlist(strsplit(as.character(x), split = ''))[1] )),
                          unlist(lapply(merged$Celltype_2, function (x) unlist(strsplit(as.character(x), split = ''))[1] )),
                          sep = '')
merged$Celltype_1_group <- unlist(lapply(merged$Celltype_1, function (x) unlist(strsplit(as.character(x), split = '\\|'))[1] ))
merged$Celltype_2_group <- unlist(lapply(merged$Celltype_2, function (x) unlist(strsplit(as.character(x), split = '\\|'))[1] ))
merged
#write.table(merged, 'heatmap_noPlt_v2.txt', quote = F, sep = '\t', row.names = F, col.names = T)


### merge ###
conervation_heatmap <- data.frame(matrix(nrow = nrow(gsva_flytofish) + nrow(gsva_flytomouse) + nrow(gsva_flytohuman), 
                                         ncol = ncol(gsva_flytofish) + ncol(gsva_fishtomouse) + ncol(gsva_mousetohuman)) )
rownames(conervation_heatmap) <- c(rownames(gsva_flytofish), rownames(gsva_flytomouse), rownames(gsva_flytohuman))
colnames(conervation_heatmap) <- c(c("Fly|PH 1", "Fly|PH", "Fly|PM", "Fly|LM", "Fly|PM 120", "Fly|Adipohemocyte", "Fly|CC", "Fly|GST-rich"), 
                                   colnames(gsva_fishtomouse), colnames(gsva_mousetohuman))

conervation_heatmap[is.na(conervation_heatmap)] <- 0
conervation_heatmap[rownames(gsva_flytofish), rownames(gsva_flytofish)] <- NA
conervation_heatmap[rownames(gsva_flytomouse), rownames(gsva_flytomouse)] <- NA
conervation_heatmap[rownames(gsva_flytofish), rownames(gsva_flytomouse)] <- NA

for (i in c(1:nrow(merged)) ) {
  ct1 <- merged[i, 'Celltype_1']
  ct2 <- merged[i, 'Celltype_2']
  conervation_heatmap[ct2, ct1] <- merged[i, 'Avgscore']
}
conervation_heatmap

conervation_heatmap_w <- data.frame(Celltype = rownames(conervation_heatmap), conervation_heatmap, check.rows = F, check.names = F)
#write.table(conervation_heatmap_w, 'heatmap_noPlt_v2.heatmap.txt', quote = F, sep = '\t', row.names = F, col.names = T)


library(pheatmap)
library(RColorBrewer)

cols <- colorRampPalette(brewer.pal(9, "OrRd"))
pheatmap(conervation_heatmap, border_color = NA, fontsize = 9, cluster_rows = F, cluster_cols = F, 
         gaps_row = c(6, 19), gaps_col = c(8, 14), na_col = 'white',
         breaks = c(0:10)*0.1, color = cols(10),
         cellheight = 8, cellwidth = 8)#,
         #filename = 'heatmap_noPlt_v2.pdf')
dev.off()


