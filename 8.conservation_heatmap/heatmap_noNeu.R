
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
# MCA
gsva_flytomouse_mca <- read.delim('../5.enrichmenttest/1.flytomouse/orthologous_allgenes/gsva_flytomouse_celltype.mca.txt', row.names = 1, check.names = F)
colnames(gsva_flytomouse_mca) <- unlist(lapply(colnames(gsva_flytomouse_mca), function (x) paste0(c('Fly|', as.character(x)), collapse = '') ))
rownames(gsva_flytomouse_mca) <- unlist(lapply(rownames(gsva_flytomouse_mca), function (x) unlist(strsplit(as.character(x), split = ' \\('))[1] ))
rownames(gsva_flytomouse_mca) <- unlist(lapply(rownames(gsva_flytomouse_mca), function (x) paste0(c('Mouse|', as.character(x)), collapse = '') ))
head(gsva_flytomouse_mca)
# TM 10X
gsva_flytomouse_tm10x <- read.delim('../5.enrichmenttest/1.flytomouse/orthologous_allgenes/gsva_flytomouse_celltype.tm10x.txt', row.names = 1, check.names = F)
colnames(gsva_flytomouse_tm10x) <- unlist(lapply(colnames(gsva_flytomouse_tm10x), function (x) paste0(c('Fly|', as.character(x)), collapse = '') ))
rownames(gsva_flytomouse_tm10x) <- unlist(lapply(rownames(gsva_flytomouse_tm10x), function (x) unlist(strsplit(as.character(x), split = ' \\('))[1] ))
rownames(gsva_flytomouse_tm10x) <- unlist(lapply(rownames(gsva_flytomouse_tm10x), function (x) paste0(c('Mouse|', as.character(x)), collapse = '') ))
head(gsva_flytomouse_tm10x)
# TM SS2
gsva_flytomouse_tmss2 <- read.delim('../5.enrichmenttest/1.flytomouse/orthologous_allgenes/gsva_flytomouse_celltype.tmss2.txt', row.names = 1, check.names = F)
colnames(gsva_flytomouse_tmss2) <- unlist(lapply(colnames(gsva_flytomouse_tmss2), function (x) paste0(c('Fly|', as.character(x)), collapse = '') ))
rownames(gsva_flytomouse_tmss2) <- unlist(lapply(rownames(gsva_flytomouse_tmss2), function (x) unlist(strsplit(as.character(x), split = ' \\('))[1] ))
rownames(gsva_flytomouse_tmss2) <- unlist(lapply(rownames(gsva_flytomouse_tmss2), function (x) paste0(c('Mouse|', as.character(x)), collapse = '') ))
head(gsva_flytomouse_tmss2)
# merged 
unique(c(rownames(gsva_flytomouse_mca), 
         unique(c(rownames(gsva_flytomouse_tm10x), 
                  rownames(gsva_flytomouse_tmss2))))) # 13

gsva_flytomouse <- data.frame(matrix(nrow = 13, ncol = ncol(gsva_flytomouse_tm10x)))
colnames(gsva_flytomouse) <- colnames(gsva_flytomouse_tm10x) # fly cells
rownames(gsva_flytomouse) <- c("Mouse|Progenitors", "Mouse|Granulocytopoietic cells", "Mouse|Erythroids", "Mouse|Neutrophils", 
                               "Mouse|Monocytes", "Mouse|Macrophage", "Mouse|DC", "Mouse|pDC", "Mouse|Basophil", 
                               "Mouse|T cells", "Mouse|NK cells", "Mouse|B cells", "Mouse|Plasma cells")

three <- intersect(rownames(gsva_flytomouse_mca), 
                   intersect(rownames(gsva_flytomouse_tm10x), rownames(gsva_flytomouse_tmss2))); three
mca_tm10x <- intersect(rownames(gsva_flytomouse_mca),
                       setdiff(rownames(gsva_flytomouse_tm10x), rownames(gsva_flytomouse_tmss2))); mca_tm10x
mca_tmss2 <- intersect(rownames(gsva_flytomouse_mca),
                       setdiff(rownames(gsva_flytomouse_tmss2), rownames(gsva_flytomouse_tm10x))); mca_tmss2
tm10x_ss2 <- intersect(rownames(gsva_flytomouse_tm10x),
                       setdiff(rownames(gsva_flytomouse_tmss2), rownames(gsva_flytomouse_mca))); tm10x_ss2
mca_unique <- setdiff( rownames(gsva_flytomouse_mca), c(rownames(gsva_flytomouse_tm10x), rownames(gsva_flytomouse_tmss2))); mca_unique

gsva_flytomouse[three, ] <- (gsva_flytomouse_mca[three, ] + gsva_flytomouse_tm10x[three, ] + gsva_flytomouse_tmss2[three, ])/3
gsva_flytomouse[mca_tm10x, ] <- (gsva_flytomouse_mca[mca_tm10x, ] + gsva_flytomouse_tm10x[mca_tm10x, ])/2
gsva_flytomouse[mca_tmss2, ] <- (gsva_flytomouse_mca[mca_tmss2, ] + gsva_flytomouse_tmss2[mca_tmss2, ])/2
gsva_flytomouse[tm10x_ss2, ] <- (gsva_flytomouse_tm10x[tm10x_ss2, ] + gsva_flytomouse_tmss2[tm10x_ss2, ])/2
gsva_flytomouse[mca_unique, ] <- gsva_flytomouse_mca[mca_unique, ]
gsva_flytomouse ###

for (ct in colnames(gsva_flytomouse)) {
  gsva_flytomouse[, ct] <- (gsva_flytomouse[, ct]-min(gsva_flytomouse[, ct]))/(max(gsva_flytomouse[, ct])-min(gsva_flytomouse[, ct]))
}
gsva_flytomouse
#################


### Fly-Human ###
# HCL
gsva_flytohuman_hcl <- read.delim('../5.enrichmenttest/1.flytohuman/orthologous_allgenes/gsva_flytohuman_celltype.hcl.txt', row.names = 1, check.names = F)
colnames(gsva_flytohuman_hcl) <- unlist(lapply(colnames(gsva_flytohuman_hcl), function (x) paste0(c('Fly|', as.character(x)), collapse = '') ))
rownames(gsva_flytohuman_hcl) <- unlist(lapply(rownames(gsva_flytohuman_hcl), function (x) unlist(strsplit(as.character(x), split = ' \\('))[1] ))
rownames(gsva_flytohuman_hcl) <- unlist(lapply(rownames(gsva_flytohuman_hcl), function (x) paste0(c('Human|', as.character(x)), collapse = '') ))
rownames(gsva_flytohuman_hcl)[11] <- 'Human|Plasma cells'
head(gsva_flytohuman_hcl)

# HCA
gsva_flytohuman_hca <- read.delim('../5.enrichmenttest/1.flytohuman/orthologous_allgenes/gsva_flytohuman_celltype.hca.txt', row.names = 1, check.names = F)
colnames(gsva_flytohuman_hca) <- unlist(lapply(colnames(gsva_flytohuman_hca), function (x) paste0(c('Fly|', as.character(x)), collapse = '') ))
rownames(gsva_flytohuman_hca) <- unlist(lapply(rownames(gsva_flytohuman_hca), function (x) unlist(strsplit(as.character(x), split = ' \\('))[1] ))
rownames(gsva_flytohuman_hca) <- unlist(lapply(rownames(gsva_flytohuman_hca), function (x) paste0(c('Human|', as.character(x)), collapse = '') ))
head(gsva_flytohuman_hca)

# merged 
unique(c( rownames(gsva_flytohuman_hcl), rownames(gsva_flytohuman_hca) )) # 17

gsva_flytohuman <- data.frame(matrix(nrow = 17, ncol = ncol(gsva_flytohuman_hcl)))
colnames(gsva_flytohuman) <- colnames(gsva_flytohuman_hcl)
rownames(gsva_flytohuman) <- c("Human|Progenitors", "Human|Granulocyte progenitor", "Human|Erythroids", "Human|Neutrophils", "Human|Monocytes", "Human|Macrophage", "Human|DC", "Human|pDC", 
                               "Human|T cells", "Human|CD8 T cells", "Human|NK cells", "Human|pre-PC", "Human|pro-B cells", "Human|pre-B cells", "Human|B cells", "Human|Plasma cells", "Human|Platelet")

both <- intersect(rownames(gsva_flytohuman_hcl), rownames(gsva_flytohuman_hca)); both
hcl_unique <- setdiff(rownames(gsva_flytohuman_hcl), rownames(gsva_flytohuman_hca)); hcl_unique
hca_unique <- setdiff(rownames(gsva_flytohuman_hca), rownames(gsva_flytohuman_hcl)); hca_unique

gsva_flytohuman[both, ] <- (gsva_flytohuman_hcl[both, ] + gsva_flytohuman_hca[both, ])/2
gsva_flytohuman[hcl_unique, ] <- gsva_flytohuman_hcl[hcl_unique, ]
gsva_flytohuman[hca_unique, ] <- gsva_flytohuman_hca[hca_unique, ]
gsva_flytohuman

for (ct in colnames(gsva_flytohuman)) {
  gsva_flytohuman[, ct] <- (gsva_flytohuman[, ct]-min(gsva_flytohuman[, ct]))/(max(gsva_flytohuman[, ct])-min(gsva_flytohuman[, ct]))
}
gsva_flytohuman
##################


### Fish-Mouse ###
# MCA
gsva_fishtomouse_mca <- read.delim('../5.enrichmenttest/2.fishtomouse/noNeu_orthologous_allgenes/gsva_fishtomouse_celltype.mca.txt', row.names = 1, check.names = F)
colnames(gsva_fishtomouse_mca) <- unlist(lapply(colnames(gsva_fishtomouse_mca), function (x) paste0(c('Fish|', as.character(x)), collapse = '') ))
rownames(gsva_fishtomouse_mca) <- unlist(lapply(rownames(gsva_fishtomouse_mca), function (x) unlist(strsplit(as.character(x), split = ' \\('))[1] ))
rownames(gsva_fishtomouse_mca) <- unlist(lapply(rownames(gsva_fishtomouse_mca), function (x) paste0(c('Mouse|', as.character(x)), collapse = '') ))
head(gsva_fishtomouse_mca)
# TM 10X
gsva_fishtomouse_tm10x <- read.delim('../5.enrichmenttest/2.fishtomouse/orthologous_allgenes/gsva_fishtomouse_celltype.tm10x.txt', row.names = 1, check.names = F)
colnames(gsva_fishtomouse_tm10x) <- unlist(lapply(colnames(gsva_fishtomouse_tm10x), function (x) paste0(c('Fish|', as.character(x)), collapse = '') ))
rownames(gsva_fishtomouse_tm10x) <- unlist(lapply(rownames(gsva_fishtomouse_tm10x), function (x) unlist(strsplit(as.character(x), split = ' \\('))[1] ))
rownames(gsva_fishtomouse_tm10x) <- unlist(lapply(rownames(gsva_fishtomouse_tm10x), function (x) paste0(c('Mouse|', as.character(x)), collapse = '') ))
head(gsva_fishtomouse_tm10x)
# TM SS2
gsva_fishtomouse_tmss2 <- read.delim('../5.enrichmenttest/2.fishtomouse/orthologous_allgenes/gsva_fishtomouse_celltype.tmss2.txt', row.names = 1, check.names = F)
colnames(gsva_fishtomouse_tmss2) <- unlist(lapply(colnames(gsva_fishtomouse_tmss2), function (x) paste0(c('Fish|', as.character(x)), collapse = '') ))
rownames(gsva_fishtomouse_tmss2) <- unlist(lapply(rownames(gsva_fishtomouse_tmss2), function (x) unlist(strsplit(as.character(x), split = ' \\('))[1] ))
rownames(gsva_fishtomouse_tmss2) <- unlist(lapply(rownames(gsva_fishtomouse_tmss2), function (x) paste0(c('Mouse|', as.character(x)), collapse = '') ))
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
##################


### Fish-Human ###
# HCL
gsva_fishtohuman_hcl <- read.delim('../5.enrichmenttest/2.fishtohuman/orthologous_allgenes/gsva_fishtohuman_celltype.hcl.txt', row.names = 1, check.names = F)
colnames(gsva_fishtohuman_hcl) <- unlist(lapply(colnames(gsva_fishtohuman_hcl), function (x) paste0(c('Fish|', as.character(x)), collapse = '') ))
rownames(gsva_fishtohuman_hcl) <- unlist(lapply(rownames(gsva_fishtohuman_hcl), function (x) unlist(strsplit(as.character(x), split = ' \\('))[1] ))
rownames(gsva_fishtohuman_hcl) <- unlist(lapply(rownames(gsva_fishtohuman_hcl), function (x) paste0(c('Human|', as.character(x)), collapse = '') ))
rownames(gsva_fishtohuman_hcl)[11] <- 'Human|Plasma cells'
head(gsva_fishtohuman_hcl)

# HCA
gsva_fishtohuman_hca <- read.delim('../5.enrichmenttest/2.fishtohuman/orthologous_allgenes/gsva_fishtohuman_celltype.hca.txt', row.names = 1, check.names = F)
colnames(gsva_fishtohuman_hca) <- unlist(lapply(colnames(gsva_fishtohuman_hca), function (x) paste0(c('Fish|', as.character(x)), collapse = '') ))
rownames(gsva_fishtohuman_hca) <- unlist(lapply(rownames(gsva_fishtohuman_hca), function (x) unlist(strsplit(as.character(x), split = ' \\('))[1] ))
rownames(gsva_fishtohuman_hca) <- unlist(lapply(rownames(gsva_fishtohuman_hca), function (x) paste0(c('Human|', as.character(x)), collapse = '') ))
head(gsva_fishtohuman_hca)

# merged 
unique(c( rownames(gsva_fishtohuman_hcl), rownames(gsva_fishtohuman_hca) )) # 17

gsva_fishtohuman <- data.frame(matrix(nrow = 17, ncol = ncol(gsva_fishtohuman_hcl)))
colnames(gsva_fishtohuman) <- colnames(gsva_fishtohuman_hcl)
rownames(gsva_fishtohuman) <- c("Human|Progenitors", "Human|Granulocyte progenitor", "Human|Erythroids", "Human|Neutrophils", "Human|Monocytes", "Human|Macrophage", "Human|DC", "Human|pDC", 
                                "Human|T cells", "Human|CD8 T cells", "Human|NK cells", "Human|pre-PC", "Human|pro-B cells", "Human|pre-B cells", "Human|B cells", "Human|Plasma cells", "Human|Platelet")

both <- intersect(rownames(gsva_fishtohuman_hcl), rownames(gsva_fishtohuman_hca)); both
hcl_unique <- setdiff(rownames(gsva_fishtohuman_hcl), rownames(gsva_fishtohuman_hca)); hcl_unique
hca_unique <- setdiff(rownames(gsva_fishtohuman_hca), rownames(gsva_fishtohuman_hcl)); hca_unique

gsva_fishtohuman[both, ] <- (gsva_fishtohuman_hcl[both, ] + gsva_fishtohuman_hca[both, ])/2
gsva_fishtohuman[hcl_unique, ] <- gsva_fishtohuman_hcl[hcl_unique, ]
gsva_fishtohuman[hca_unique, ] <- gsva_fishtohuman_hca[hca_unique, ]
gsva_fishtohuman

for (ct in colnames(gsva_fishtohuman)) {
  gsva_fishtohuman[, ct] <- (gsva_fishtohuman[, ct]-min(gsva_fishtohuman[, ct]))/(max(gsva_fishtohuman[, ct])-min(gsva_fishtohuman[, ct]))
}
gsva_fishtohuman
###################


### Mouse-Human ###
# HCL
gsva_mousetohuman_hcl <- read.delim('../5.enrichmenttest/3.mousetohuman/noNeu_orthologous_allgenes/gsva_mousetohuman_celltype.hcl.txt', row.names = 1, check.names = F)
colnames(gsva_mousetohuman_hcl) <- unlist(lapply(colnames(gsva_mousetohuman_hcl), function (x) paste0(c('Mouse|', as.character(x)), collapse = '') ))
rownames(gsva_mousetohuman_hcl) <- unlist(lapply(rownames(gsva_mousetohuman_hcl), function (x) unlist(strsplit(as.character(x), split = ' \\('))[1] ))
rownames(gsva_mousetohuman_hcl) <- unlist(lapply(rownames(gsva_mousetohuman_hcl), function (x) paste0(c('Human|', as.character(x)), collapse = '') ))
rownames(gsva_mousetohuman_hcl)[11] <- 'Human|Plasma cells'
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


meta_flytomouse <- read.delim('../6.MetaNeighbor/1.flytomouse/celltype/MetaNeighbor_flytomouse_orthologous.celltype_top_hits.txt')
for (i in c(1:nrow(meta_flytomouse)) ) {
  ct1 <- as.character(meta_flytomouse[i, 1])
  ct2 <- as.character(meta_flytomouse[i, 2])
  gsvascore <- c(gsva_flytomouse[ct1, ct2], gsva_flytomouse[ct2, ct1])
  meta_flytomouse[i, 'avgscore'] <- mean(c(meta_flytomouse[i,3], gsvascore))
}
head(meta_flytomouse)


meta_flytohuman <- read.delim('../6.MetaNeighbor/1.flytohuman/celltype/MetaNeighbor_flytohuman_orthologous.celltype_top_hits.txt')
for (i in c(1:nrow(meta_flytohuman)) ) {
  ct1 <- as.character(meta_flytohuman[i, 1])
  ct2 <- as.character(meta_flytohuman[i, 2])
  gsvascore <- c(gsva_flytohuman[ct1, ct2], gsva_flytohuman[ct2, ct1])
  meta_flytohuman[i, 'avgscore'] <- mean(c(meta_flytohuman[i,3], gsvascore))
}
head(meta_flytohuman)


### Fish ###
meta_fishtomouse <- read.delim('../6.MetaNeighbor/2.fishtomouse_noNeu/celltype_merged/MetaNeighbor_fishtomouse_orthologous.celltype_top_hits.txt')
for (i in c(1:nrow(meta_fishtomouse)) ) {
  ct1 <- as.character(meta_fishtomouse[i, 1])
  ct2 <- as.character(meta_fishtomouse[i, 2])
  gsvascore <- c(gsva_fishtomouse[ct1, ct2], gsva_fishtomouse[ct2, ct1])
  meta_fishtomouse[i, 'avgscore'] <- mean(c(meta_fishtomouse[i,3], gsvascore))
}
head(meta_fishtomouse)

meta_fishtohuman <- read.delim('../6.MetaNeighbor/2.fishtohuman/celltype_merged/MetaNeighbor_fishtohuman_orthologous.celltype_top_hits.txt')
for (i in c(1:nrow(meta_fishtohuman)) ) {
  ct1 <- as.character(meta_fishtohuman[i, 1])
  ct2 <- as.character(meta_fishtohuman[i, 2])
  gsvascore <- c(gsva_fishtohuman[ct1, ct2], gsva_fishtohuman[ct2, ct1])
  meta_fishtohuman[i, 'avgscore'] <- mean(c(meta_fishtohuman[i,3], gsvascore))
}
head(meta_fishtohuman)


### Mouse ###
meta_mousetohuman <- read.delim('../6.MetaNeighbor/3.mousetohuman_noNeu/celltype_merged/MetaNeighbor_orthologous.celltype_top_hits.txt')
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
merged[c(7, 14), c(1,2)] <- merged[c(7, 14), c(2,1)]
merged[c(18, 20, 23), c(1,2)] <- merged[c(18, 20, 23), c(2,1)]
merged[c(29, 30, 32, 36, 38:40, 42), c(1,2)] <- merged[c(29, 30, 32, 36, 38:40, 42), c(2,1)]
merged[c(48, 51, 54, 55, 57:60), c(1,2)] <- merged[c(48, 51, 54, 55, 57:60), c(2,1)]
merged

merged$hierarchy <- paste(unlist(lapply(merged$Celltype_1, function (x) unlist(strsplit(as.character(x), split = ''))[1] )),
                          unlist(lapply(merged$Celltype_2, function (x) unlist(strsplit(as.character(x), split = ''))[1] )),
                          sep = '')
merged$Celltype_1_group <- unlist(lapply(merged$Celltype_1, function (x) unlist(strsplit(as.character(x), split = '\\|'))[1] ))
merged$Celltype_2_group <- unlist(lapply(merged$Celltype_2, function (x) unlist(strsplit(as.character(x), split = '\\|'))[1] ))
merged
#write.table(merged, 'heatmap_noNeu.txt', quote = F, sep = '\t', row.names = F, col.names = T)


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
#write.table(conervation_heatmap_w, 'heatmap_noNeu.heatmap.txt', quote = F, sep = '\t', row.names = F, col.names = T)


library(pheatmap)
library(RColorBrewer)

cols <- colorRampPalette(brewer.pal(9, "OrRd"))
pheatmap(conervation_heatmap, border_color = NA, fontsize = 9, cluster_rows = F, cluster_cols = F, 
         gaps_row = c(6, 19), gaps_col = c(8, 14), na_col = 'white',
         breaks = c(0:10)*0.1, color = cols(10),
         cellheight = 8, cellwidth = 8)#,
         #filename = 'heatmap_noNeu.pdf')
dev.off()




