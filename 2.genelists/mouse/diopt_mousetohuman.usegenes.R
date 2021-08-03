
mousetohuman <- read.delim('diopt_mousetohuman.txt', check.names = F)
head(mousetohuman)

### human expression union set ###
exprs_human <- readRDS('../../../Model_species/human/merged_commongenes/tmp/expr.pseudobulk.Rds')
nrow(exprs_human) # 20935
exprs_human$sums <- rowSums(exprs_human)
exprs_human <- subset(exprs_human, sums != 0)
exprs_human$sums <- NULL
nrow(exprs_human) # 20131

genes_human <- rownames(exprs_human); length(genes_human) # 20131 genes


### mouse expression union set ###
exprs_mouse <- readRDS('../../../Model_species/mouse/merged_commongenes_v3/tmp/mouse_pseudobulk.Rds')
head(exprs_mouse)
nrow(exprs_mouse) # 11504
exprs_mouse$sums <- rowSums(exprs_mouse)
exprs_mouse <- subset(exprs_mouse, sums != 0)
exprs_mouse$sums <- NULL
nrow(exprs_mouse) # 11504

genes_mouse <- rownames(exprs_mouse); length(genes_mouse) # 11504 genes


###
mousetohuman_best <- subset(mousetohuman, BestScore == 'Yes')
mousetohuman_best <- subset(mousetohuman_best, DIOPTScore > 1)
mousetohuman_best <- subset(mousetohuman_best, HumanSymbol %in% genes_human)
mousetohuman_best <- subset(mousetohuman_best, MouseGeneID %in% genes_mouse)
mousetohuman_best <- droplevels(mousetohuman_best)
nrow(mousetohuman_best); length(unique(mousetohuman_best$MouseGeneID)) # 10646 unique genes in mouse

freq <- data.frame(table(mousetohuman_best$MouseGeneID))
freq <- subset(freq, Freq > 1); head(freq)

singleton_best <- subset(mousetohuman_best, !MouseGeneID %in% freq$Var1); nrow(singleton_best) # 10515

### one-to-many ###
selected_best <- data.frame(matrix(ncol = ncol(singleton_best), nrow = 0))
colnames(selected_best) <- colnames(singleton_best)
evenscores <- selected_best
undefined <- c()

for (gene in freq$Var1) {
  tmp <- subset(mousetohuman_best, MouseGeneID == gene)
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

nrow(singleton_best) # 10515
nrow(selected_best) # 96
nrow(evenscores) # 0

### many-to-one ###
use_best <- rbind(singleton_best, selected_best)
nrow(use_best); length(unique(use_best$HumanSymbol)) # 10611 mouse genes mapped to 10380 human genes

freq <- data.frame(table(use_best$HumanSymbol))
freq <- subset(freq, Freq > 1); head(freq)

onetoone <- subset(use_best, !HumanSymbol %in% freq$Var1); nrow(onetoone) # 10234

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
    tmp_exprs <- exprs_mouse[as.character(tmp$MouseGeneID), ]
    tmp_exprs$Sums <- rowSums(tmp_exprs)
    tmp_exprs_hi <- rownames(tmp_exprs[tmp_exprs$Sums == max(tmp_exprs$Sums), ])
    print (tmp_exprs_hi)
    
    if (length(tmp_exprs_hi) == 1) {
      tmp <- subset(tmp, MouseGeneID == tmp_exprs_hi)
      manytoone_selected <- rbind(manytoone_selected, tmp)
    } else {
      tmp <- subset(tmp, MouseGeneID %in% tmp_exprs_hi)
      evenscores <- rbind(evenscores, tmp)
    }
  }
}
nrow(onetoone) # 10234
nrow(manytoone_selected) # 146

orthologous_genes <- rbind(onetoone, manytoone_selected)
nrow(orthologous_genes) # 10380
length(unique(orthologous_genes$MouseGeneID)) # 10380
length(unique(orthologous_genes$HumanSymbol)) # 10380

#write.table(orthologous_genes, 'diopt_mousetohuman.usegenes.txt', row.names = F, col.names = T, quote = F, sep = '\t')


