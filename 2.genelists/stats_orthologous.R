
flytofish <- read.delim('fly/diopt_flytofish.usegenes.txt', check.names = F)
head(flytofish); nrow(flytofish) # 5676
flytomouse <- read.delim('fly/diopt_flytomouse.usegenes.txt', check.names = F)
head(flytomouse); nrow(flytomouse) # 5119
flytohuman <- read.delim('fly/diopt_flytohuman.usegenes.txt', check.names = F)
head(flytohuman); nrow(flytohuman) # 6385

fishtomouse <- read.delim('fish/diopt_fishtomouse.usegenes.txt', check.names = F)
head(fishtomouse); nrow(fishtomouse) # 8713

mousetohuman <- read.delim('mouse/diopt_mousetohuman.usegenes.txt', check.names = F)
head(mousetohuman); nrow(mousetohuman) # 10380


flytofishtomouse <- intersect(flytofish$ZebrafishSymbol, fishtomouse$ZebrafishGeneID); length(flytofishtomouse) # 4266
length(intersect(fishtomouse[fishtomouse$ZebrafishGeneID %in% flytofishtomouse, 'MouseSymbol'], # == 4266
                 flytomouse$MouseSymbol)) # 3805 genes intersected

flytofishtomousetohuman <- intersect(fishtomouse[fishtomouse$ZebrafishGeneID %in% intersect(flytofish$ZebrafishSymbol, fishtomouse$ZebrafishGeneID) , 'MouseSymbol'],
                                     mousetohuman$MouseGeneID); length(flytofishtomousetohuman) # 4219

mouse_flytofishtomousetohuman <- mousetohuman[mousetohuman$MouseGeneID %in% flytofishtomousetohuman,]
head(mouse_flytofishtomousetohuman); nrow(mouse_flytofishtomousetohuman)
fish_flytofishtomousetohuman <- fishtomouse[fishtomouse$MouseSymbol %in% mouse_flytofishtomousetohuman$MouseGeneID, ]
head(fish_flytofishtomousetohuman); nrow(fish_flytofishtomousetohuman)
fly_flytofishtomousetohuman <- flytofish[flytofish$ZebrafishSymbol %in% fish_flytofishtomousetohuman$ZebrafishGeneID, ]
head(fly_flytofishtomousetohuman); nrow(fly_flytofishtomousetohuman)

merged <- fly_flytofishtomousetohuman[, c(1,8)]

rownames(merged) <- merged$ZebrafishSymbol
rownames(fish_flytofishtomousetohuman) <- as.character(fish_flytofishtomousetohuman$ZebrafishGeneID)
fish_flytofishtomousetohuman <- fish_flytofishtomousetohuman[as.character(merged$ZebrafishSymbol), ]
merged <- cbind(merged, fish_flytofishtomousetohuman[, c(1,8)])

rownames(merged) <- merged$MouseSymbol
rownames(mouse_flytofishtomousetohuman) <- as.character(mouse_flytofishtomousetohuman$MouseGeneID)
mouse_flytofishtomousetohuman <- mouse_flytofishtomousetohuman[as.character(merged$MouseSymbol), ]
merged <- cbind(merged, mouse_flytofishtomousetohuman[, c(1,8)])

merged$ZebrafishSymbol <- NULL
merged$MouseSymbol <- NULL
head(merged); nrow(merged)

#write.table(merged, 'stats_orthologous.4219genes.txt', quote = F, sep = '\t', row.names = F, col.names = T)


length(intersect(mousetohuman[mousetohuman$MouseGeneID %in% flytofishtomousetohuman, 'HumanSymbol'], # == 4219
                 flytohuman$HumanSymbol)) # 3766 genes intersected

