library(wordcloud)
library(dplyr)

#replace sentence in a row
pathway_counts$Pathway <- gsub('HSFA7/ HSFA6B-regulatory network-induced by drought and ABA. ', 
                             'Regulatory network-induced by drought and ABA. ', pathway_counts$Pathway)

# Data processing: Count the number of unique genes per pathway
pathway_counts <- ann_top_blast %>%
  group_by(Pathway) %>%
  summarise(GeneCount = n_distinct(OsID))

# Creating the word cloud
color_palette <- c("blue", "olivedrab")

png("wordcloud.png", width = 100, height = 200)
wordcloud(words = sum_cgene$chr, freq = sum_cgene$count_gene_per_chr, 
          min.freq = 1, random.order = FALSE, rot.per = 0.1, 
          colors = color_palette)
dev.off()

sum_cgene$chr <- factor(sum_cgene$chr, levels = sort(unique(sum_cgene$chr))) #sort value

ggplot(sum_cgene, aes(x=chr, y=count_gene_per_chr)) +
  geom_segment(aes(x=chr, xend=chr, y=0, yend=count_gene_per_chr), color="skyblue") +
  geom_point(color="blue", size=4, alpha=0.6) +
  theme_light() +
  coord_flip() +
  theme(
    panel.grid.major.y = element_blank(),
    panel.border = element_blank(),
    axis.ticks.y = element_blank()
  )



***************************

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

#cloudword
library(wordcloud2)
wordcloud2(pvalue_fltr, color = "random-light", backgroundColor = "grey")

###data <- data %>%
  rownames_to_column(var="the name you want")
