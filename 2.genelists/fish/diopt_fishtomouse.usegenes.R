
fishtomouse <- read.delim('diopt_fishtomouse.txt', check.names = F)
head(fishtomouse)

### mouse expression union set ###
exprs_mouse <- readRDS('../../../Model_species/mouse/merged_commongenes_v3/tmp/mouse_pseudobulk.Rds')
head(exprs_mouse)
nrow(exprs_mouse) # 11504
exprs_mouse$sums <- rowSums(exprs_mouse)
exprs_mouse <- subset(exprs_mouse, sums != 0)
exprs_mouse$sums <- NULL
nrow(exprs_mouse) # 11504

genes_mouse <- rownames(exprs_mouse); length(genes_mouse) # 11504 genes


exprs_fish <- readRDS('../../../Model_species/zebrafish/Tang_Q_et_al/InDrops/tmp/celltype.pseudocell10.Rds')
exprs_fish$sums <- rowSums(exprs_fish)
exprs_fish <- subset(exprs_fish, sums != 0)
exprs_fish$sums <- NULL
genes_fish <- rownames(exprs_fish); length(genes_fish) # 22522 genes


fishtomouse_best <- subset(fishtomouse, BestScore == 'Yes')
fishtomouse_best <- subset(fishtomouse_best, DIOPTScore > 1)
fishtomouse_best <- subset(fishtomouse_best, MouseSymbol %in% genes_mouse)
fishtomouse_best <- subset(fishtomouse_best, ZebrafishGeneID %in% genes_fish)
fishtomouse_best <- droplevels(fishtomouse_best)
nrow(fishtomouse_best); length(unique(fishtomouse_best$ZebrafishGeneID)) # 11020 unique genes in fly

freq <- data.frame(table(fishtomouse_best$ZebrafishGeneID))
freq <- subset(freq, Freq > 1); head(freq)

singleton_best <- subset(fishtomouse_best, !ZebrafishGeneID %in% freq$Var1); nrow(singleton_best) # 10730

### one-to-many ###
selected_best <- data.frame(matrix(ncol = ncol(singleton_best), nrow = 0))
colnames(selected_best) <- colnames(singleton_best)
evenscores <- selected_best
undefined <- c()

for (gene in freq$Var1) {
  tmp <- subset(fishtomouse_best, ZebrafishGeneID == gene)
  maxscore <- max(tmp$WeightedScore)
  tmp <- subset(tmp, WeightedScore == maxscore)
  if (nrow(tmp) != 1){
    tmp <- subset(tmp, MouseSymbol != '')
    if (nrow(tmp) == 1) {
      selected_best <- rbind(selected_best, tmp_2)
    } else {
      tmp_2 <- subset(tmp, BestScoreReverse == 'Yes')
      if (nrow(tmp_2) == 1) {
        selected_best <- rbind(selected_best, tmp_2)
      } else {
        ### consider expression ###
        tmp_exprs <- exprs_mouse[as.character(tmp_2$MouseSymbol), ]
        tmp_exprs$Sums <- rowSums(tmp_exprs)
        tmp_exprs_hi <- rownames(tmp_exprs[tmp_exprs$Sums == max(tmp_exprs$Sums), ])
        print (tmp_exprs_hi)
        
        if (length(tmp_exprs_hi) == 1) {
          tmp_2 <- subset(tmp_2, MouseSymbol == tmp_exprs_hi)
          selected_best <- rbind(selected_best, tmp_2)
        } else {
          tmp_2 <- subset(tmp_2, MouseSymbol %in% tmp_exprs_hi)
          evenscores <- rbind(evenscores, tmp_2)
        }
      }
    }
  } else {
    selected_best <- rbind(selected_best, tmp)
    }
  }

nrow(singleton_best) # 10730
nrow(selected_best) # 198

### many-to-one ###
use_best <- rbind(singleton_best, selected_best)
nrow(use_best); length(unique(use_best$MouseSymbol)) # 10928 fish genes mapped to 8713 mouse genes

freq <- data.frame(table(use_best$MouseSymbol))
freq <- subset(freq, Freq > 1); head(freq)

onetoone <- subset(use_best, !MouseSymbol %in% freq$Var1); nrow(onetoone) # 6763

manytoone_selected <- data.frame(matrix(ncol = ncol(use_best), nrow = 0))
colnames(manytoone_selected) <- colnames(use_best)
evenscores <- manytoone_selected
undefined <- c()

for (gene in freq$Var1) {
  
  tmp <- subset(use_best, MouseSymbol == gene)
  maxscore <- max(tmp$WeightedScore)
  tmp <- subset(tmp, WeightedScore == maxscore)
  
  if (nrow(tmp) == 1) {
    manytoone_selected <- rbind(manytoone_selected, tmp)
  } else {
    
    ### consider expression ###
    tmp_exprs <- exprs_fish[as.character(tmp$ZebrafishGeneID), ]
    tmp_exprs$Sums <- rowSums(tmp_exprs)
    tmp_exprs_hi <- rownames(tmp_exprs[tmp_exprs$Sums == max(tmp_exprs$Sums), ])
    print (tmp_exprs_hi)
    
    if (length(tmp_exprs_hi) == 1) {
      tmp <- subset(tmp, ZebrafishGeneID == tmp_exprs_hi)
      manytoone_selected <- rbind(manytoone_selected, tmp)
    } else {
      tmp <- subset(tmp, ZebrafishGeneID %in% tmp_exprs_hi)
      evenscores <- rbind(evenscores, tmp)
    }
  }
}
nrow(onetoone) # 6763
nrow(manytoone_selected) # 1950

orthologous_genes <- rbind(onetoone, manytoone_selected)
nrow(orthologous_genes) # 8713
length(unique(orthologous_genes$ZebrafishGeneID)) # 8713
length(unique(orthologous_genes$MouseSymbol)) # 8713

#write.table(orthologous_genes, 'diopt_fishtomouse.usegenes.txt', row.names = F, col.names = T, quote = F, sep = '\t')


