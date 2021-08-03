
fishtohuman <- read.delim('diopt_fishtohuman.txt', check.names = F)
head(fishtohuman)

### human expression union set ###
exprs_human <- readRDS('../../../Model_species/human/merged_commongenes/tmp/expr.pseudobulk.Rds')
nrow(exprs_human) # 20935
exprs_human$sums <- rowSums(exprs_human)
exprs_human <- subset(exprs_human, sums != 0)
exprs_human$sums <- NULL
nrow(exprs_human) # 20131
genes_human <- rownames(exprs_human); length(genes_human) # 20131 genes


exprs_fish <- readRDS('../../../Model_species/zebrafish/Tang_Q_et_al/InDrops/tmp/celltype.pseudocell10.Rds')
exprs_fish$sums <- rowSums(exprs_fish)
exprs_fish <- subset(exprs_fish, sums != 0)
exprs_fish$sums <- NULL
genes_fish <- rownames(exprs_fish); length(genes_fish) # 22522 genes


fishtohuman_best <- subset(fishtohuman, BestScore == 'Yes')
fishtohuman_best <- subset(fishtohuman_best, DIOPTScore > 1)
fishtohuman_best <- subset(fishtohuman_best, HumanSymbol %in% genes_human)
fishtohuman_best <- subset(fishtohuman_best, ZebrafishGeneID %in% genes_fish)
fishtohuman_best <- droplevels(fishtohuman_best)
nrow(fishtohuman_best); length(unique(fishtohuman_best$ZebrafishGeneID)) # 15307 unique genes in fly

freq <- data.frame(table(fishtohuman_best$ZebrafishGeneID))
freq <- subset(freq, Freq > 1); head(freq)

singleton_best <- subset(fishtohuman_best, !ZebrafishGeneID %in% freq$Var1); nrow(singleton_best) # 14822

### one-to-many ###
selected_best <- data.frame(matrix(ncol = ncol(singleton_best), nrow = 0))
colnames(selected_best) <- colnames(singleton_best)
evenscores <- selected_best
undefined <- c()

for (gene in freq$Var1) {
  tmp <- subset(fishtohuman_best, ZebrafishGeneID == gene)
  maxscore <- max(tmp$WeightedScore)
  tmp <- subset(tmp, WeightedScore == maxscore)
  if (nrow(tmp) != 1){
    tmp <- subset(tmp, HumanSymbol != '')
    if (nrow(tmp) == 1) {
      selected_best <- rbind(selected_best, tmp_2)
    } else {
      tmp_2 <- subset(tmp, BestScoreReverse == 'Yes')
      if (nrow(tmp_2) == 1) {
        selected_best <- rbind(selected_best, tmp_2)
      } else {
        ### consider expression ###
        tmp_exprs <- exprs_human[as.character(tmp_2$HumanSymbol), ]
        tmp_exprs$Sums <- rowSums(tmp_exprs)
        tmp_exprs_hi <- rownames(tmp_exprs[tmp_exprs$Sums == max(tmp_exprs$Sums), ])
        print (tmp_exprs_hi)
        
        if (length(tmp_exprs_hi) == 1) {
          tmp_2 <- subset(tmp_2, HumanSymbol == tmp_exprs_hi)
          selected_best <- rbind(selected_best, tmp_2)
        } else {
          tmp_2 <- subset(tmp_2, HumanSymbol %in% tmp_exprs_hi)
          evenscores <- rbind(evenscores, tmp_2)
        }
      }
    }
  } else {
    selected_best <- rbind(selected_best, tmp)
    }
  }

nrow(singleton_best) # 14822
nrow(selected_best) # 342
nrow(evenscores) # 

### many-to-one ###
use_best <- rbind(singleton_best, selected_best)
nrow(use_best); length(unique(use_best$HumanSymbol)) # 15164 fish genes mapped to 12044 human genes

freq <- data.frame(table(use_best$HumanSymbol))
freq <- subset(freq, Freq > 1); head(freq)

onetoone <- subset(use_best, !HumanSymbol %in% freq$Var1); nrow(onetoone) # 9285

manytoone_selected <- data.frame(matrix(ncol = ncol(use_best), nrow = 0))
colnames(manytoone_selected) <- colnames(use_best)
evenscores <- manytoone_selected
undefined <- c()

for (gene in freq$Var1) {
  
  tmp <- subset(use_best, HumanSymbol == gene)
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
nrow(onetoone) # 9285
nrow(manytoone_selected) # 2759

orthologous_genes <- rbind(onetoone, manytoone_selected)
nrow(orthologous_genes) # 12044
length(unique(orthologous_genes$ZebrafishGeneID)) # 12044
length(unique(orthologous_genes$HumanSymbol)) # 12044

#write.table(orthologous_genes, 'diopt_fishtohuman.usegenes.txt', row.names = F, col.names = T, quote = F, sep = '\t')


