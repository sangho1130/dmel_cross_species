
flytomouse <- read.delim('diopt_flytomouse.txt', check.names = F)
head(flytomouse)

### mouse expression union set ###
exprs_mouse <- readRDS('../../../Model_species/mouse/merged_commongenes_v3/tmp/mouse_pseudobulk.Rds')
nrow(exprs_mouse) # 11504
exprs_mouse$sums <- rowSums(exprs_mouse)
exprs_mouse <- subset(exprs_mouse, sums != 0)
exprs_mouse$sums <- NULL
nrow(exprs_mouse) # 11504
genes_mouse <- rownames(exprs_mouse); length(genes_mouse) # 11504 genes


exprs_fly <- readRDS('../../../Model_species/fly/tmp/celltype_v2_pseudocell10.Rds')
exprs_fly$sums <- rowSums(exprs_fly)
exprs_fly <- subset(exprs_fly, sums != 0)
exprs_fly$sums <- NULL
genes_fly <- rownames(exprs_fly); length(genes_fly) # 16175 genes


flytomouse_best <- subset(flytomouse, BestScore == 'Yes')
flytomouse_best <- subset(flytomouse_best, DIOPTScore > 1)
flytomouse_best <- subset(flytomouse_best, MouseSymbol %in% genes_mouse)
flytomouse_best <- subset(flytomouse_best, FlyGeneID %in% genes_fly)
flytomouse_best <- droplevels(flytomouse_best)
nrow(flytomouse_best); length(unique(flytomouse_best$FlyGeneID)) # 6852 unique genes in fly

freq <- data.frame(table(flytomouse_best$FlyGeneID))
freq <- subset(freq, Freq > 1); head(freq)

singleton_best <- subset(flytomouse_best, !FlyGeneID %in% freq$Var1); nrow(singleton_best) # 6089

### one-to-many ###
selected_best <- data.frame(matrix(ncol = ncol(singleton_best), nrow = 0))
colnames(selected_best) <- colnames(singleton_best)
evenscores <- selected_best
undefined <- c()

for (gene in freq$Var1) {
  tmp <- subset(flytomouse_best, FlyGeneID == gene)
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

nrow(singleton_best) # 6089
nrow(selected_best) # 454

### many-to-one ###
use_best <- rbind(singleton_best, selected_best)
nrow(use_best); length(unique(use_best$MouseSymbol)) # 6543 fly genes mapped to 5119 mouse genes

freq <- data.frame(table(use_best$MouseSymbol))
freq <- subset(freq, Freq > 1); head(freq)

onetoone <- subset(use_best, !MouseSymbol %in% freq$Var1); nrow(onetoone) # 4406

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
nrow(onetoone) # 4406
nrow(manytoone_selected) # 713

orthologous_genes <- rbind(onetoone, manytoone_selected)
nrow(orthologous_genes) # 5119
length(unique(orthologous_genes$FlyGeneID)) # 5119
length(unique(orthologous_genes$MouseSymbol)) # 5119

#write.table(orthologous_genes, 'diopt_flytomouse.usegenes.txt', row.names = F, col.names = T, quote = F, sep = '\t')


