#making graph
plot(sscore, type= "b", lwd= 0.8, col = "blue",
     xlab = "inflation", ylab = "SScore",
     main = "SScore")

p2 <- ggplot(sscore, aes(inflation, SScore)) + geom_line() + geom_point()
p2 + scale_y_continuous("SScore",breaks = c(2,4,6,8,10,12,14)) +
  scale_x_continuous("inflation", 
                     breaks = c(1.0, 1.1, 1.2,
                                1.3, 1.4, 1.5,
                                1.6, 1.7, 1.8,
                                1.9, 2.0, 2.1,
                                2.2, 2.3, 2.4,
                                2.5, 2.6, 2.7, 
                                2.8, 2.9, 3.0))
                                       

colnames(i12)[2]="ID"
nic_i12<- anti_join(i12, go_i12, by="ID")

#group a & b
i15_clstr_hits<- merge(go_i15, hits_terms, by="ID")
i15_clstr_nohits<- anti_join(go_i15, i15_clstr_hits, by="ID")

#group c & d
i15_nic_hits<- merge(nic_i15, nic_hits_terms, by="ID")#c
i15_nic_nohit<- anti_join(nic_i15, i15_nic_hits, by="ID")#d

#library(stats)
df <- data.frame("a" = c(735, 278), "b" = c(262, 27), 
                 row.names = c("c", "d"))

fisher.test(df)

#analysis
hub<- read.csv("~/Desktop/mcl_i22_node.csv")
gg<- read.csv("~/Desktop/mcl_i22_edge.csv")

top_hub_ann<- top_hub_ann %>%
  filter(ID %in% c("Em_1059",
"Em_3169",
"Em_1266",
"Em_2844",
"gene2030",
"gene6362",
"gene7661",
"gene1579",
"gene5057",
"gene2000"))

b<- a %>% pivot_longer(cols= -LocID, names_to="SRA", values_to="FPKM")
