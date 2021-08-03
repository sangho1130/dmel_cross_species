
flytofish <- read.delim('diopt_flytofish.txt', check.names = F)
head(flytofish)

exprs_fish <- readRDS('../../../Model_species/zebrafish/Tang_Q_et_al/InDrops/tmp/celltype.pseudocell10.Rds'); head(exprs_fish)
exprs_fish$sums <- rowSums(exprs_fish)
exprs_fish <- subset(exprs_fish, sums != 0)
exprs_fish$sums <- NULL
genes_fish <- rownames(exprs_fish); length(genes_fish) # 22522 genes

exprs_fly <- readRDS('../../../Model_species/fly/tmp/celltype_v2_pseudocell10.Rds')
exprs_fly$sums <- rowSums(exprs_fly)
exprs_fly <- subset(exprs_fly, sums != 0)
exprs_fly$sums <- NULL
genes_fly <- rownames(exprs_fly); length(genes_fly) # 16175 genes


flytofish_best <- subset(flytofish, BestScore == 'Yes')
flytofish_best <- subset(flytofish_best, DIOPTScore > 1)
flytofish_best <- subset(flytofish_best, ZebrafishSymbol %in% genes_fish)
flytofish_best <- subset(flytofish_best, FlyGeneID %in% genes_fly)
flytofish_best <- droplevels(flytofish_best)
nrow(flytofish_best); length(unique(flytofish_best$FlyGeneID)) # 7811 unique genes in fly

freq <- data.frame(table(flytofish_best$FlyGeneID))
freq <- subset(freq, Freq > 1); head(freq)

singleton_best <- subset(flytofish_best, !FlyGeneID %in% freq$Var1); nrow(singleton_best) # 6368

### one-to-many ###
selected_best <- data.frame(matrix(ncol = ncol(singleton_best), nrow = 0))
colnames(selected_best) <- colnames(singleton_best)
evenscores <- selected_best
undefined <- c()

for (gene in freq$Var1) {
  tmp <- subset(flytofish_best, FlyGeneID == gene)
  maxscore <- max(tmp$WeightedScore)
  tmp <- subset(tmp, WeightedScore == maxscore)
  if (nrow(tmp) != 1){
    tmp <- subset(tmp, ZFINID != '')
    if (nrow(tmp) == 1) {
      selected_best <- rbind(selected_best, tmp_2)
    } else {
      tmp_2 <- subset(tmp, BestScoreReverse == 'Yes')
      if (nrow(tmp_2) == 1) {
        selected_best <- rbind(selected_best, tmp_2)
      } else {
        ### consider expression ###
        tmp_exprs <- exprs_fish[as.character(tmp_2$ZebrafishSymbol), ]
        tmp_exprs$Sums <- rowSums(tmp_exprs)
        tmp_exprs_hi <- rownames(tmp_exprs[tmp_exprs$Sums == max(tmp_exprs$Sums), ])
        print (tmp_exprs_hi)
        
        if (length(tmp_exprs_hi) == 1) {
          tmp_2 <- subset(tmp_2, ZebrafishSymbol == tmp_exprs_hi)
          selected_best <- rbind(selected_best, tmp_2)
        } else {
          tmp_2 <- subset(tmp_2, ZebrafishSymbol %in% tmp_exprs_hi)
          evenscores <- rbind(evenscores, tmp_2)
        }
      }
    }
  } else {
    selected_best <- rbind(selected_best, tmp)
  }
}

nrow(singleton_best) # 6368
nrow(selected_best) # 752
nrow(evenscores)

### many-to-one ###
use_best <- rbind(singleton_best, selected_best)
nrow(use_best); length(unique(use_best$ZebrafishSymbol)) # 7120 fly genes mapped to 5676 fish genes

freq <- data.frame(table(use_best$ZebrafishSymbol))
freq <- subset(freq, Freq > 1); head(freq)

onetoone <- subset(use_best, !ZebrafishSymbol %in% freq$Var1); nrow(onetoone) # 4940

manytoone_selected <- data.frame(matrix(ncol = ncol(use_best), nrow = 0))
colnames(manytoone_selected) <- colnames(use_best)
evenscores <- manytoone_selected
undefined <- c()

for (gene in freq$Var1) {
  
  tmp <- subset(use_best, ZebrafishSymbol == gene)
  maxscore <- max(tmp$WeightedScore)
  tmp <- subset(tmp, WeightedScore == maxscore)
  
  if (nrow(tmp) == 1) {
    manytoone_selected <- rbind(manytoone_selected, tmp)
  } else {
    
    ### consider expression ###
    tmp_exprs <- exprs_fly[as.character(tmp$FlyGeneID), ]
    tmp_exprs$Sums <- rowSums(tmp_exprs)
    tmp_exprs_hi <- rownames(tmp_exprs[tmp_exprs$Sums == max(tmp_exprs$Sums), ])
    print (tmp_exprs_hi)
    
    if (length(tmp_exprs_hi) == 1) {
      tmp <- subset(tmp, FlyGeneID == tmp_exprs_hi)
      manytoone_selected <- rbind(manytoone_selected, tmp)
    } else {
      tmp <- subset(tmp, FlyGeneID %in% tmp_exprs_hi)
      evenscores <- rbind(evenscores, tmp)
    }
  }
}
nrow(onetoone) # 4940
nrow(manytoone_selected) # 736

orthologous_genes <- rbind(onetoone, manytoone_selected)
nrow(orthologous_genes) # 5676
length(unique(orthologous_genes$FlyGeneID)) # 5676
length(unique(orthologous_genes$ZebrafishSymbol)) # 5676

#write.table(orthologous_genes, 'diopt_flytofish.usegenes.txt', row.names = F, col.names = T, quote = F, sep = '\t')


