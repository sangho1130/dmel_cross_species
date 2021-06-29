
mousetohuman <- read.delim('diopt_mousetohuman.txt', check.names = F)
head(mousetohuman)

### human expression union set ###
hcl <- readRDS('../../human/hcl/tmp/celltype_refined.pseudocell10.Rds')
colnames(hcl) <- unlist( lapply(colnames(hcl), function (x) paste0(c(as.character(x), '(HCL)'), collapse = ' ') ) )
head(hcl); dim(hcl)

hca <- readRDS('../../human/hca/tmp/expr.pseudobulk.Rds')
colnames(hca) <- unlist( lapply(colnames(hca), function (x) paste0(c(as.character(x), '(HCA)'), collapse = ' ') ) )
head(hca); dim(hca)

union <- unique(c(rownames(hcl), rownames(hca)))
exprs_human <- data.frame(matrix(nrow = length(union), ncol = ncol(hcl) + ncol(hca)))
colnames(exprs_human) <- c(colnames(hcl), colnames(hca)); rownames(exprs_human) <- union
exprs_human[is.na(exprs_human)] <- 0

exprs_human[rownames(hcl), colnames(hcl)] <- hcl
exprs_human[rownames(hca), colnames(hca)] <- hca
head(exprs_human)
exprs_human$sums <- rowSums(exprs_human)
exprs_human <- subset(exprs_human, sums != 0)
exprs_human$sums <- NULL
genes_human <- rownames(exprs_human); length(genes_human) # 45089 genes


### mouse expression union set ###
mca <- readRDS('../../mouse/mca/tmp/celltype_refined.pseudocell10.Rds')
mca <- readRDS('../../mouse/mca/tmp/celltype_refined_noNeu.pseudocell10.Rds')
colnames(mca) <- unlist( lapply(colnames(mca), function (x) paste0(c(as.character(x), '(MCA)'), collapse = ' ') ) )
head(mca); dim(mca)

tabulamuris_10x <- readRDS('../../mouse/tabulamuris_10x/tmp/celltype_refined.pseudocell10.Rds')
colnames(tabulamuris_10x) <- unlist( lapply(colnames(tabulamuris_10x), function (x) paste0(c(as.character(x), '(TM10X)'), collapse = ' ') ) )
head(tabulamuris_10x); dim(tabulamuris_10x)

tabulamuris_ss2 <- readRDS('../../mouse/tabulamuris_ss2/tmp/celltype_refined.pseudocell10.Rds')
colnames(tabulamuris_ss2) <- unlist( lapply(colnames(tabulamuris_ss2), function (x) paste0(c(as.character(x), '(TMSS2)'), collapse = ' ') ) )
head(tabulamuris_ss2); dim(tabulamuris_ss2)


union <- unique(c(rownames(mca), rownames(tabulamuris_10x), rownames(tabulamuris_ss2)))
exprs_mouse <- data.frame(matrix( nrow = length(union), ncol = ncol(mca) + ncol(tabulamuris_10x) + ncol(tabulamuris_ss2) ))
colnames(exprs_mouse) <- c(colnames(mca), colnames(tabulamuris_10x), colnames(tabulamuris_ss2)); rownames(exprs_mouse) <- union
exprs_mouse[is.na(exprs_mouse)] <- 0

exprs_mouse[rownames(mca), colnames(mca)] <- mca
exprs_mouse[rownames(tabulamuris_10x), colnames(tabulamuris_10x)] <- tabulamuris_10x
exprs_mouse[rownames(tabulamuris_ss2), colnames(tabulamuris_ss2)] <- tabulamuris_ss2
head(exprs_mouse)
genes_mouse$sums <- rowSums(genes_mouse)
genes_mouse <- subset(genes_mouse, sums != 0)
genes_mouse$sums <- NULL

genes_mouse <- rownames(exprs_mouse); length(genes_mouse) # 24455 genes


mousetohuman_best <- subset(mousetohuman, BestScore == 'Yes')
mousetohuman_best <- subset(mousetohuman_best, DIOPTScore > 1)
mousetohuman_best <- subset(mousetohuman_best, HumanSymbol %in% genes_human)
mousetohuman_best <- subset(mousetohuman_best, MouseGeneID %in% genes_mouse)
mousetohuman_best <- droplevels(mousetohuman_best)
nrow(mousetohuman_best); length(unique(mousetohuman_best$MouseGeneID)) # 17968 unique genes in mouse

freq <- data.frame(table(mousetohuman_best$MouseGeneID))
freq <- subset(freq, Freq > 1); head(freq)

singleton_best <- subset(mousetohuman_best, !MouseGeneID %in% freq$Var1); nrow(singleton_best) # 17551

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

nrow(singleton_best) # 17551
nrow(selected_best) # 248
nrow(evenscores) # 

### many-to-one ###
use_best <- rbind(singleton_best, selected_best)
nrow(use_best); length(unique(use_best$HumanSymbol)) # 17799 mouse genes mapped to 15840 human genes

freq <- data.frame(table(use_best$HumanSymbol))
freq <- subset(freq, Freq > 1); head(freq)

onetoone <- subset(use_best, !HumanSymbol %in% freq$Var1); nrow(onetoone) # 14403

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
nrow(onetoone) # 14403
nrow(manytoone_selected) # 1437

orthologous_genes <- rbind(onetoone, manytoone_selected)
nrow(orthologous_genes) # 15840
length(unique(orthologous_genes$MouseGeneID)) # 15840
length(unique(orthologous_genes$HumanSymbol)) # 15840

#write.table(orthologous_genes, 'diopt_mousetohuman.usegenes.txt', row.names = F, col.names = T, quote = F, sep = '\t')


