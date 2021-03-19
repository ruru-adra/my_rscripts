library(pheatmap)
library(tibble)

a<- a %>% column_to_rownames("geneID")
a<- as.matrix(a)

#z-score function
cal_z_score <- function(x){
  (x - mean(x)) / sd(x)
}

#calculate z_score
a2 <- t(apply(a, 1, cal_z_score))
pheatmap(a)

