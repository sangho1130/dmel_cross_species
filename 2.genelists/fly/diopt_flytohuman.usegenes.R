
flytohuman <- read.delim('diopt_flytohuman.txt', check.names = F)
head(flytohuman)

### mouse expression union set ###
exprs_human <- readRDS('../../../Model_species/human/merged_commongenes/tmp/expr.pseudobulk.Rds')
nrow(exprs_human) # 20935
exprs_human$sums <- rowSums(exprs_human)
exprs_human <- subset(exprs_human, sums != 0)
exprs_human$sums <- NULL
nrow(exprs_human) # 20131
genes_human <- rownames(exprs_human); length(genes_human) # 20131 genes


exprs_fly <- readRDS('../../../Model_species/fly/tmp/celltype_v2_pseudocell10.Rds')
exprs_fly$sums <- rowSums(exprs_fly)
exprs_fly <- subset(exprs_fly, sums != 0)
exprs_fly$sums <- NULL
genes_fly <- rownames(exprs_fly); length(genes_fly) # 16175 genes


flytohuman_best <- subset(flytohuman, BestScore == 'Yes')
flytohuman_best <- subset(flytohuman_best, DIOPTScore > 1)
flytohuman_best <- subset(flytohuman_best, HumanSymbol %in% genes_human)
flytohuman_best <- subset(flytohuman_best, FlyGeneID %in% genes_fly)
flytohuman_best <- droplevels(flytohuman_best)
nrow(flytohuman_best); length(unique(flytohuman_best$FlyGeneID)) # 8650 unique genes in fly

freq <- data.frame(table(flytohuman_best$FlyGeneID))
freq <- subset(freq, Freq > 1); head(freq)

singleton_best <- subset(flytohuman_best, !FlyGeneID %in% freq$Var1); nrow(singleton_best) # 7271

### one-to-many ###
selected_best <- data.frame(matrix(ncol = ncol(singleton_best), nrow = 0))
colnames(selected_best) <- colnames(singleton_best)
evenscores <- selected_best
undefined <- c()

for (gene in freq$Var1) {
  tmp <- subset(flytohuman_best, FlyGeneID == gene)
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

nrow(singleton_best) # 7271
nrow(selected_best) # 820
nrow(evenscores) # 

### many-to-one ###
use_best <- rbind(singleton_best, selected_best)
nrow(use_best); length(unique(use_best$HumanSymbol)) # 8091 fly genes mapped to 6385 mouse genes

freq <- data.frame(table(use_best$HumanSymbol))
freq <- subset(freq, Freq > 1); head(freq)

onetoone <- subset(use_best, !HumanSymbol %in% freq$Var1); nrow(onetoone) # 5464

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
nrow(onetoone) # 5464
nrow(manytoone_selected) # 921

orthologous_genes <- rbind(onetoone, manytoone_selected)
nrow(orthologous_genes) # 6385
length(unique(orthologous_genes$FlyGeneID)) # 6385
length(unique(orthologous_genes$HumanSymbol)) # 6385

#write.table(orthologous_genes, 'diopt_flytohuman.usegenes.txt', row.names = F, col.names = T, quote = F, sep = '\t')


