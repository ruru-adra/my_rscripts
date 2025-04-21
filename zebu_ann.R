library(dplyr)
library(tidyr)
library(tidyverse)
library(stringi)
library(stringr)
library(tibble)
library(data.table)
library(sqldf)

num_snp<- clean_sort_genic %>%
  filter(!chr %in% c("X", "Y")) %>%
  arrange(as.numeric(chr), pos)

ysnp<- clean_sort_genic %>%
  filter(chr=="Y") %>%
  arrange(pos)

sort_genic<- rbind(num_snp, x_snp, y_snp)

sum_genic<- sqldf("select count(snpID) as count_snp, snpID
                  from sort_genic group by snpID
                  order by count_snp DESC")

snp_more_than_one_gene<- merge(snp_more_than_one_gene,
                               sort_genic,by="snpID")

clean_snp_more_than_one_gene<- snp_more_than_one_gene %>%
  group_by(snpID) %>%
  filter(n() > 1) %>%
  slice(1) %>%
  ungroup()

snp_gene_sgle<- sum_genic %>%
  filter(count_snp==1) %>%
  merge(sort_genic, by="snpID")

clean_sort_genic<- rbind(clean_snp_more_than_one_gene, snp_gene_sgle)

num_snp_ann<- num_snp %>%
  merge(mrna_zebu, by="tID") %>%
  select(snpID, chr, pos, ref, alel, tID, protID, gene_desc, gene_symbol, locus_tag,
         feature_type, eff, cDNA_pos,
         CDS_pos, prot_pos, aa, codons)

xsnp_ann<- xsnp %>%
  merge(mrna_zebu, by="tID") %>%
  select(snpID, chr, pos, ref, alel, tID, protID, gene_desc, gene_symbol, locus_tag,
         feature_type, eff, cDNA_pos,
         CDS_pos, prot_pos, aa, codons)

ysnp_ann<- ysnp %>%
  merge(mrna_zebu, by="tID") %>%
  select(snpID, chr, pos, ref, alel, tID, protID, gene_desc, gene_symbol, locus_tag,
         feature_type, eff, cDNA_pos,
         CDS_pos, prot_pos, aa, codons)
         
num_snp_ann<- num_snp_ann %>%
           arrange(as.numeric(chr), pos)

ysnp_ann<- ysnp_ann %>%
  arrange(pos)
         
clean_sort_genic_ann<- rbind(num_snp_ann, xsnp_ann, ysnp_ann)         


write.table(clean_sort_genic_ann, "clean_sort_genic_ann.txt", col.names=T, row.names=F, sep="\t", quote=F)
write.table(ysnp_ann, "ysnp_ann.txt", col.names=T, row.names=F, sep="\t", quote=F)

#total annotated snp #12,622,674 snp
sum_genic<- sqldf("select count(snpID) as count_snp, snpID
                  from clean_sort_genic_ann group by snpID
                  order by count_snp DESC")

#annotated snp by chr
library("CMplot")

CMplot(clean_sort_genic_ann,plot.type="d",bin.size=1e10,
       chr.den.col=c("darkgreen", "yellow", "red"),
       file="jpg",file.name="",dpi=300,
       main="zebu wgs",file.output=TRUE,verbose=TRUE,width=9,height=6)

#annotated genic snp categories
snp_data<- sqldf("select count(eff) as Total_SNPs, chr, eff from clean_sort_genic_ann
          group by chr, eff")

sum_eff<- sqldf("select count(eff) as Total_SNPs, eff from clean_sort_genic_ann
          group by eff")


library(pheatmap)
snp_matrix <- tidyr::pivot_wider(snp_data, names_from = eff, values_from = Total_SNPs, values_fill = 0)
snp_matrix <- as.data.frame(snp_matrix)
rownames(snp_matrix) <- snp_matrix$chr
snp_matrix <- snp_matrix[,-1]

pheatmap(as.matrix(snp_matrix), color = colorRampPalette(c("darkgreen", "yellow", "red"))(100),
         cluster_rows = FALSE, cluster_cols = FALSE, main = "SNP Effects Across Chromosomes")

#
chr<- read.table("chr_zebu.txt", header=T, sep="\t")
tsnp<- read.table("tsnp_zebu_final.txt", header=T, sep="\t")

colnames(tsnp)[2]="refseq_chr"

new_tsnp<- merge(chr, tsnp, by="refseq_chr")
  
new_tsnp<- select(snpID, chr_a, chr, pos, allele, NELORE, BRAHMAN, BRAKMAS, CHOLASTANI, HAINAN, KK1, KK2, KKE, LEIQONG, MIRKADIM,
         RCC, RS, SAHIWAL)

new_tsnp<- paste(new_tsnp$chr, new_tsnp$pos, sep="_")

num_snp<- new_tsnp %>%
  filter(!chr %in% c("X", "Y"))

xsnp<- new_tsnp %>%
  filter(chr %in% c("X"))

ysnp<- new_tsnp %>%
  filter(chr %in% c("Y"))
  
a_select<- select(a, snpID, chr, pos, allele, NELORE, BRAHMAN, BRAKMAS, CHOLASTANI, HAINAN, KK1, KK2, KKE, LEIQONG, MIRKADIM, RCC, RS, SAHIWAL)

#chr plot for tsnp
library(CMplot)
CMplot(for_map, plot.type="d", bin.size=1e6,
       chr.den.col=c("darkgreen", "yellow", "red"),
       file="jpg",file.name="",dpi=300,
       main="zebu wgs",file.output=TRUE,verbose=TRUE,width=9,height=6)

#total intergenic vs genic
#38,371,918 total snp
#25,749,244 intergenic/not annotated into genic


#
setDT(chr)
uq_snp<- read.table("unique_snp.txt", header=T, sep="\t")
setDT(uq_snp)

uq_snp<- uq_snp[, .(a,b,pos)]
colnames(uq_snp)[2]="split_snp"

uq_snp[, c("a","b","pos") := tstrsplit(snpID, "_")]

uq_snp<- uq_snp[, refseq_chr := paste(a, b, sep = "_")] #paste dt

uq_snp_merge<- merge(chr, uq_snp, by="refseq_chr")
uq_snp_merge<- uq_snp_merge[, snpID := paste(chr_a, pos, sep="_")]

uq_snp_merge<- uq_snp_merge[, .(snpID)]
uq_snp_merge<- as.data.frame(uq_snp_merge)

setDT(uq_snp_merge)
a2<- merge(uq_snp_merge,a, by="snpID")


#compare brakmaks, KK1, KK2 and KKE
snp_local<- select(tsnp, snpID, NELORE, BRAKMAS, KK1, KK2, KKE)
colnames(snp_local)[2]="ref"

#select snp that identical to ref #12612671
snp_identical<- snp_local %>%
  filter(if_all(c(BRAKMAS, KK1, KK2, KKE), ~ . == ref))
  
#total number snp brakmas, kk1, kk2, kke $25,759,247
tsnp_local<- anti_join(snp_local, snp_identical, by="snpID")

#tsnp brakmas
bms_tsnp<- sqldf("select snpID, ref, BRAKMAS from tsnp_local
                where BRAKMAS != REF")#14,795,632

kke_tsnp<- sqldf("select snpID, ref, KKE from tsnp_local
                where KKE != REF")#15,167,294

kk1_tsnp<- sqldf("select snpID, ref, KK1 from tsnp_local
                where KK1 != REF")#10756767

kk2_tsnp<- sqldf("select snpID, ref, KK2 from tsnp_local
                where KK2 != REF")#10900372

#brakmas=3,873,811
bms_uq_snp<- sqldf("select snpID, ref, BRAKMAS from tsnp_local
                where BRAKMAS != REF
                AND BRAKMAS != KK1
                AND BRAKMAS != KK2
                AND BRAKMAS != KKE")
               
#kke=3,632,541
KKE_uq_snp<- sqldf("select snpID, ref, KKE from tsnp_local
                where KKE != REF
                AND KKE != KK1
                AND KKE != KK2
                AND KKE != BRAKMAS")


#kk1=
KK1_uq_snp<- sqldf("select snpID, ref, KK1 from tsnp_local
                where KK1 != REF
                AND KK1 != KKE
                AND KK1 != KK2
                AND KK1 != BRAKMAS")


#kk2=
KK2_uq_snp<- sqldf("select snpID, ref, KK2 from tsnp_local
                where KK2 != REF
                AND KK2 != KK1
                AND KK2 != KKE
                AND KK2 != BRAKMAS")




#homo & hetero


#transition and transversion ratio
is_tr<- with(kke_tsnp, (ref == "A" & KKE == "G") | 
               (ref == "G" & KKE  == "A") |
               (ref == "C" & KKE  == "T") | 
               (ref == "T" & KKE  == "C"))

is_tv <- with(kke_tsnp, 
                        (ref == "A" & KKE %in% c("C", "T")) |
                          (ref == "G" & KKE %in% c("C", "T")) |
                          (ref == "C" & KKE %in% c("A", "G")) |
                          (ref == "T" & KKE %in% c("A", "G")))

# Count transitions and transversions
tr <- sum(is_tr, na.rm = TRUE)
tv <- sum(is_tv, na.rm = TRUE)

Ts_Tv_ratio <- tr/tv


#
is_ht <- with(kk2hthm, ref != KK2)

# Identify homozygous SNPs (Ref == Alt)
is_hm<- with(kk2hthm, ref == KK2)

# Count the SNPs
ht_count <- sum(is_ht, na.rm = TRUE)
hm_count <- sum(is_hm, na.rm = TRUE)

#ann
kk1_ann<- merge(kk1_tsnp, ann_snp, by="snpID")
sum_kk1_ann<- sqldf("select count(snpID) as count_snp, eff from kk1_ann
                    group by eff")

kk2_ann<- merge(kk2_tsnp, ann_snp, by="snpID")
sum_kk2_ann<- sqldf("select count(snpID) as count_snp, eff from kk2_ann
                    group by eff")


#
kk1_breed<- select(kk1_ann, snpID, gene_symbol)
kk1_breed$breed="kk1"


unique_genes_per_breed <- eff_breed %>%
  group_by(gene_symbol) %>% 
  filter(n_distinct(breed) == 1) %>%  # Keep genes appearing in only one breed
  ungroup()

# Split results by breed
unique_genes_list <- split(unique_genes_per_breed, unique_genes_per_breed$breed)

# Save results for each breed
for (breed in names(unique_genes_list)) {
  write.csv(unique_genes_list[[breed]], paste0("Unique_Genes_", breed, ".csv"), row.names = FALSE)
}

#deleterious
del_snpID<- select(del_snp, snpID)
del_bms<-  merge(bms_ann, del_snpID, by="snpID")

#
to_binary <- to_binary %>%
  mutate(binary = ifelse(ref == alel, 0, 1))

# Convert to wide format (Binary Matrix)
breed_binary <-to_binary %>%
  select(-ref, -alel) %>%
  pivot_wider(names_from = breed, values_from = binary, values_fill = list(binary = 0))


##
library(ggvenn)
library(ggplot2)
breed_list <- list(
  Brakmas = rownames(breed_binary)[breed_binary$Brakmas == 1],
  KKE = rownames(breed_binary)[breed_binary$KKE == 1],
  KK1 = rownames(breed_binary)[breed_binary$KK1 == 1],
  KK2 = rownames(breed_binary)[breed_binary$KK2 == 1]
)

venn_plot <- ggvenn(
  breed_list, 
  fill_color = c("#E63946", "#457B9D", "#2A9D8F", "#F4A261")
)

# Save the plot
ggsave("venn_diagram.png", venn_plot, width = 6, height = 6, dpi = 300)  


#
bms_eco_snp<- sqldf("select snpID, ref, alel, BRAKMAS, KK1, KK2, KKE,
trait, gene_symbol, eff, chr, pos from breed_snp_eco
                where BRAKMAS != REF
                AND BRAKMAS != KK1
                AND BRAKMAS != KK2
                AND BRAKMAS != KKE")

kke_eco_snp<- sqldf("select snpID, ref, alel, BRAKMAS, KK1, KK2, KKE,
trait, gene_symbol, eff, chr, pos from breed_snp_eco
                where KKE != REF
                AND KKE != KK1
                AND KKE != KK2
                AND KKE != BRAKMAS")


kk1_eco_snp<- sqldf("select snpID, ref, alel, BRAKMAS, KK1, KK2, KKE,
trait, gene_symbol, eff, chr, pos from breed_snp_eco
                where KK1 != REF
                AND KK1 != KKE
                AND KK1 != KK2
                AND KK1 != BRAKMAS")


breed_eco_uq_snp<- rbind(bms_eco_snp, KK1_uq_snp, kk2_eco_snp, kke_eco_snp)



#######plot del snp
ggplot(fltr_del_to_plot, aes(x = total_snp, y = reorder(gene_symbol, total_snp),
                             fill = gene_symbol)) +
  geom_bar(stat = "identity") +  # Bar plot
  theme_minimal() +  # Clean theme
  labs(
    title = "Top 10 Genes with Most Deleterious SNPs",
    x = "Total Number of SNPs",
    y = "Gene Symbol"
  ) +
  scale_fill_brewer(palette = "Set3") +  # Color palette
  theme(axis.text.x = element_text(size = 12), 
        axis.text.y = element_text(size = 12))  # Adjust text size


#