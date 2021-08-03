library(Seurat)
library(plyr)
library(dplyr)
library(ggplot2)

myObj <- readRDS('../../Drop-seq_merge_harmony_6.22/tmp/dropseq.combined.harmony_flt.Rds')
head(myObj@meta.data)

summary(myObj@meta.data$subclustering)
myObj <- subset(myObj, celltype != 'PSC')
myObj@meta.data <- droplevels(myObj@meta.data)

### Cell type v1 ###
myObj@meta.data$celltype_v1 <- mapvalues(myObj@meta.data$subclustering,
                                         from = levels(myObj@meta.data$subclustering),
                                         to = c("PH 1", "PH", "PH", "PH", "PH", "PH", "PM", "PM 120", "PM 120", "PM 120", 
                                                "LM", "LM", "CC", "CC", "GST-rich", "Adipohemocyte"))
Idents(myObj) <- 'celltype_v1'
levels(myObj@meta.data$celltype_v1)

degmarkers <- FindAllMarkers(object = myObj, min.pct = 0.25, logfc.threshold = 0.25, only.pos = T, test.use = 'MAST')
write.table(degmarkers, 'celltype_v1.txt', sep = '\t', quote = F, col.names = T, row.names = F) 
saveRDS(degmarkers, 'celltype_v1.Rds')

top10 <- degmarkers %>% group_by(cluster) %>% top_n(10, avg_logFC)
dehm <- DoHeatmap(object = myObj, features = top10$gene, angle = 90, size = 3, raster = T, draw.lines = F)
ggsave('celltype_v1.pdf', units = 'cm', width = 40, height = 30)


### Cell type v1 - Lymph Gland ###
myObj@meta.data$celltype_v1 <- mapvalues(myObj@meta.data$subclustering,
                                         from = levels(myObj@meta.data$subclustering),
                                         to = c("PH 1", "PH", "PH", "PH", "PH", "PH", "PM", "PM 120", "PM 120", "PM 120", 
                                                "LM", "LM", "CC", "CC", "GST-rich", "Adipohemocyte"))
myObj_l <- subset(myObj, origin == 'Lymph gland')
myObj_l@meta.data <- droplevels(myObj_l@meta.data)
Idents(myObj_l) <- 'celltype_v1'
levels(myObj_l@meta.data$celltype_v1)

degmarkers <- FindAllMarkers(object = myObj_l, min.pct = 0.25, logfc.threshold = 0.25, only.pos = T, test.use = 'MAST')
write.table(degmarkers, 'celltype_v1_lg.txt', sep = '\t', quote = F, col.names = T, row.names = F) 
saveRDS(degmarkers, 'celltype_v1_lg.Rds')

top10 <- degmarkers %>% group_by(cluster) %>% top_n(10, avg_logFC)
dehm <- DoHeatmap(object = myObj_l, features = top10$gene, angle = 90, size = 3, raster = T, draw.lines = F)
ggsave('celltype_v1_lg.pdf', units = 'cm', width = 30, height = 20)

### Cell type v1 - Circulation ###
myObj@meta.data$celltype_v1 <- mapvalues(myObj@meta.data$subclustering,
                                         from = levels(myObj@meta.data$subclustering),
                                         to = c("PH 1", "PH", "PH", "PH", "PH", "PH", "PM", "PM 120", "PM 120", "PM 120", 
                                                "LM", "LM", "CC", "CC", "GST-rich", "Adipohemocyte"))
myObj_c <- subset(myObj, origin != 'Lymph gland')
myObj_c@meta.data <- droplevels(myObj_c@meta.data)
Idents(myObj_c) <- 'celltype_v1'
levels(myObj_c@meta.data$celltype_v1)

degmarkers <- FindAllMarkers(object = myObj_c, min.pct = 0.25, logfc.threshold = 0.25, only.pos = T, test.use = 'MAST')
write.table(degmarkers, 'celltype_v1_circ.txt', sep = '\t', quote = F, col.names = T, row.names = F) 
saveRDS(degmarkers, 'celltype_v1_circ.Rds')

top10 <- degmarkers %>% group_by(cluster) %>% top_n(10, avg_logFC)
dehm <- DoHeatmap(object = myObj_c, features = top10$gene, angle = 90, size = 3, raster = T, draw.lines = F)
ggsave('celltype_v1_circ.pdf', units = 'cm', width = 30, height = 20)
###


