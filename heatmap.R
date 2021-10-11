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

#ggplot (https://www.r-graph-gallery.com/index.html)

#plotting GO (dot plot)
library(ggplot2)

ggplot(pvalue_fltr) +
  geom_point(aes(x= p.value, y= description, color= description))

#plotting GO (barplot)
ggplot(topgo, aes(x=description, y=p.value)) + 
  geom_bar(stat = "identity") +
  coord_flip()

#plotting barplot
pvalue_fltr %>%
  mutate(name = fct_reorder(description, p.value)) %>%
  ggplot( aes(x=description, y=p.value)) +
  geom_bar(stat="identity", fill="#f68060", alpha=.6, width=.4) +
  coord_flip() +
  xlab("") +
  theme_bw()
