mcl dr_pcc2mcl.txt -i 1.1 --abc -o dr_out11.txt
mcl ripe_pcc2mcl.txt -i 1.1 --abc -o ripe_out11.txt


#split by char pos and extract
mp$a<- substr(mp$new, 9, nchar(mp$new))


#PASTE MULTIPLE ROWS IN SINGLE ROWS
a<- main_matches %>% 
    group_by(winner) %>% 
    summarise(years = paste(year, collapse=", "))


library(tidyr)
library(dplyr)
library(stringi)
library(stringr)
library(sqldf)
library(CMplot)
library(readxl)

A0025<- read.table("CC_SNP/SNP_cc_A0025.txt", header=T, sep="\t", quote="")
B0025<- read.table("CC_SNP/SNP_cc_B0025.txt", header=T, sep="\t", quote="")
Bali<- read.table("CC_SNP/SNP_cc_Bali.txt", header=T, sep="\t", quote="")
MR12H<- read.table("CC_SNP/SNP_cc_KADARIA.txt", header=T, sep="\t", quote="")
MASRIA<- read.table("CC_SNP/SNP_cc_MASRIA.txt", header=T, sep="\t", quote="")
MH63<- read.table("CC_SNP/SNP_cc_MH63.txt", header=T, sep="\t", quote="")
MR103<- read.table("CC_SNP/SNP_cc_MR103.txt", header=T, sep="\t", quote="")
MR106<- read.table("CC_SNP/SNP_cc_MR106.txt", header=T, sep="\t", quote="")
MR127<- read.table("CC_SNP/SNP_cc_MR127.txt", header=T, sep="\t", quote="")
MR297<- read.table("CC_SNP/SNP_cc_MR297.txt", header=T, sep="\t", quote="")
MR307<- read.table("CC_SNP/SNP_cc_MR307.txt", header=T, sep="\t", quote="")
MRM16<- read.table("CC_SNP/SNP_cc_MRM16.txt", header=T, sep="\t", quote="")
MRQ100<- read.table("CC_SNP/SNP_cc_MRQ100.txt", header=T, sep="\t", quote="")
MRQ76<- read.table("CC_SNP/SNP_cc_MRQ76.txt", header=T, sep="\t", quote="")
PH9<- read.table("CC_SNP/SNP_cc_PH9.txt", header=T, sep="\t", quote="")
PM<- read.table("CC_SNP/SNP_cc_PM.txt", header=T, sep="\t", quote="")
PS<- read.table("CC_SNP/SNP_cc_PS.txt", header=T, sep="\t", quote="")
MRQ74<- read.table("CC_SNP/SNP_cc_Q74.txt", header=T, sep="\t", quote="")
R004<- read.table("CC_SNP/SNP_cc_R004.txt", header=T, sep="\t", quote="")
V20B<- read.table("CC_SNP/SNP_cc_V20B.txt", header=T, sep="\t", quote="")
ZS97<- read.table("CC_SNP/SNP_cc_ZS97.txt", header=T, sep="\t", quote="")

names(MASRIA)<- c("chr", "ref", "alel", "qual", "pos", "GT", "AD", "DP", "GQ", "PL")

fltr_BALI<- BALI %>% 
  filter(DP >= 3) %>%
  filter(GQ >= 10)

fltr_MASRIA$snpID<- paste(fltr_MASRIA$chr, fltr_MASRIA$pos, sep="_")

fltr_ZS97<- filter(fltr_ZS97, alel %in% c("A", "C", "G", "T"))

#genotype format

g_MR127<- unique(select(fltr_MR12H, snpID, ref, alel))
colnames(g_MR127)[3]="MR127"

tsnp<- Reduce(function(x,y) merge(x, y, by = "snpID", all.x = TRUE, all.y = TRUE),
               list(PM, PH9, MASRIA, MRM16, MRQ76, PS, MR297, MR307, MR12H, 
                    MR127, MRQ74, MR103, MR106)) #total snp

gsnp<- Reduce(function(x,y) merge(x, y, by = c("snpID", "ref", "alel"), all.x = TRUE, all.y = TRUE),
              list(PM, PH9, MASRIA, MRM16, MRQ76, PS, MR297, MR307, MR12H, 
                   MR127, MRQ74, MR103, MR106)) #get py allele

sum_tsnp<- sqldf("select count(snpID) AS count_alel, snpID from gsnp group by snpID")

a<- tsnp
a[is.na(a)]<- "NC"
a$countNC<- rowSums(a == "NC")


sgle<- filter(sum_tsnp, count_alel == 1)
py<- filter(sum_tsnp, count_alel != 1)

b<- py_alel %>% #pivot to get py allele in single row
  pivot_wider(names_from = alel, values_from=alel)
b$alel<- paste(b$G, b$A, b$T, b$C, sep="/")


#annotation
genic<- filter(fltr_ann, snp_effect %in% c("missense_variant", "synonymous_variant",
                                           "3_prime_UTR_variant", "5_prime_UTR_variant",
                                           "intron_variant"))
ssrg_genic  <- genic %>%
merge(ssrg, by="R498_ID") %>%
merge(hapmap, by=c("snpID", "chr", "pos", "ref", "alel")) %>%
select(snpID, chr, pos, ref, alel, snp_effect, R498_ID, category, geneName, 
       geneDesc, AA, AA_POS, A0025, B0025, BALI, MASRIA, MH63,  
       MR103, MR106, MR127, MR12H, MR297, MR307, MRM16, MRQ100,
       MRQ74, MRQ76, PH9, PM, PS, R004, V20B, ZS97, type) %>%
filter(type != "monomorphic") %>%
unique()

ssrg_del<- ssrg_genic %>%
  merge(del_sift4g, by=c("R498_ID", "snpID", "chr", "pos", "ref", "alel")) %>%
  select(snpID, chr, pos, ref, alel, snp_effect, R498_ID, TRANSCRIPT_ID, 
         category, geneName, 
         geneDesc, AA, AA_POS,  REF_AMINO, ALT_AMINO, AMINO_POS,      
         SIFT_SCORE, SIFT_MEDIAN, NUM_SEQS, SIFT_PREDICTION, A0025, B0025, BALI, 
         MASRIA, MH63,  
         MR103, MR106, MR127, MR12H, MR297, MR307, MRM16, MRQ100,
         MRQ74, MRQ76, PH9, PM, PS, R004, V20B, ZS97, type) %>%
  filter(snp_effect == "missense_variant") %>%
  unique()

rf_genic <- genic %>%
  merge(rf, by="R498_ID") %>%
  merge(hapmap, by=c("snpID", "chr", "pos", "ref", "alel")) %>%
  select(snpID, chr, pos, ref, alel, snp_effect, R498_ID, geneName, 
         geneDesc, domain, AA, AA_POS, A0025, B0025, BALI, MASRIA, MH63,  
         MR103, MR106, MR127, MR12H, MR297, MR307, MRM16, MRQ100,
         MRQ74, MRQ76, PH9, PM, PS, R004, V20B, ZS97, type) %>%
  filter(type != "monomorphic") %>%
  unique()

rf_del<- rf_genic %>%
  merge(del_sift4g, by=c("R498_ID", "snpID", "chr", "pos", "ref", "alel")) %>%
  select(snpID, chr, pos, ref, alel, snp_effect, R498_ID, geneName, 
         geneDesc, domain, AA, AA_POS, A0025, B0025, BALI, MASRIA, MH63,  
         MR103, MR106, MR127, MR12H, MR297, MR307, MRM16, MRQ100,
         MRQ74, MRQ76, PH9, PM, PS, R004, V20B, ZS97, type) %>%
  filter(snp_effect == "missense_variant") %>%
  unique()


mbkbase_genic <- genic %>%
  merge(mbkbase, by="R498_ID") %>%
  merge(hapmap, by=c("snpID", "chr", "pos", "ref", "alel")) %>%
  select(snpID, chr, pos, ref, alel, snp_effect, R498_ID, geneName, 
         geneDesc, Trait, Class, AA, AA_POS, A0025, B0025, BALI, MASRIA, MH63,  
         MR103, MR106, MR127, MR12H, MR297, MR307, MRM16, MRQ100,
         MRQ74, MRQ76, PH9, PM, PS, R004, V20B, ZS97, type) %>%
  unique()

mbkbase_del<- mbkbase_genic %>%
  merge(del_sift4g, by=c("R498_ID", "snpID", "chr", "pos", "ref", "alel")) %>%
  select(snpID, chr, pos, ref, alel, snp_effect, R498_ID, geneName, 
         geneDesc, Trait, Class, AA, AA_POS, A0025, B0025, BALI, MASRIA, MH63,  
         MR103, MR106, MR127, MR12H, MR297, MR307, MRM16, MRQ100,
         MRQ74, MRQ76, PH9, PM, PS, R004, V20B, ZS97, type) %>%
  filter(snp_effect == "missense_variant") %>%
  unique()

hyb_genic <- genic %>%
  merge(hyb, by="R498_ID") %>%
  merge(hapmap, by=c("snpID", "chr", "pos", "ref", "alel")) %>%
  select(snpID, chr, pos, ref, alel, snp_effect, R498_ID, geneName, 
         geneDesc, Traits, LocID, AA, AA_POS, A0025, B0025, MR12H,R004, type) %>%
  filter(type != "monomorphic") %>%
  unique()

hyb_del<- hyb_genic %>%
  merge(del_sift4g, by=c("R498_ID", "snpID", "chr", "pos", "ref", "alel")) %>%
  select(snpID, chr, pos, ref, alel, snp_effect, R498_ID, geneName, 
         geneDesc, Traits, LocID, AA, AA_POS, A0025, B0025, MR12H,R004, type) %>%
  filter(snp_effect == "missense_variant") %>%
  unique()
  
#pattern analysis

hyb_uq<- sqldf("select * from hyb_genic where MR12H != A0025
                    AND MR12H != B0025
                    AND MR12H != R004 ")

parent_hyb<- sqldf("select * from hyb_genic where 
                    MR12H != A0025
                    AND MR12H != B0025
                    AND MR12H != R004
                    AND A0025 == B0025
                    AND A0025 != R004
                    AND B0025 != R004 ")
                     
parent_mbkbase<- sqldf("select * from mbkbase_genic where 
                    MR12H != A0025
                    AND MR12H != B0025
                    AND MR12H != R004
                    AND A0025 == B0025
                    AND A0025 != R004
                    AND B0025 != R004 ")

parent_rf<- sqldf("select * from rf_genic where 
                    MR12H != A0025
                    AND MR12H != B0025
                    AND MR12H != R004
                    AND A0025 == B0025
                    AND A0025 != R004
                    AND B0025 != R004 ")   



waxy_pt<- sqldf("select * from hapmap where PH9 == MASRIA
                    AND PH9 == PM
                    AND MASRIA == PM
                    AND PH9 != MRM16
                    AND PH9 != MRQ76
                    AND PH9 != PS
                    AND PH9 != MR297
                    AND PH9 != MR307
                    AND PH9 != MR12H
                    AND PH9 != MR127
                    AND PH9 != MRQ74
                    AND PH9 != MR103
                    AND PH9 != MR106
                    AND MASRIA != MRM16
                    AND MASRIA != MRQ76
                    AND MASRIA != PS
                    AND MASRIA != MR297
                    AND MASRIA != MR307
                    AND MASRIA != MR12H
                    AND MASRIA != MR127
                    AND MASRIA != MRQ74
                    AND MASRIA != MR103
                    AND MASRIA != MR106 
                    AND PM != MRM16
                    AND PM != MRQ76
                    AND PM != PS
                    AND PM != MR297
                    AND PM != MR307
                    AND PM != MR12H
                    AND PM != MR127
                    AND PM != MRQ74
                    AND PM != MR103
                    AND PM != MR106 ")

low_pt<- sqldf("select * from hapmap where MRM16 == MRQ76
                    AND MRM16 == PS
                    AND MRQ76 == PS
                    AND MRM16 != PH9
                    AND MRM16 != MASRIA
                    AND MRM16 != MR297
                    AND MRM16 != MR307
                    AND MRM16 != MR12H
                    AND MRM16 != MR127
                    AND MRM16 != MRQ74
                    AND MRM16 != MR103
                    AND MRM16 != MR106
                    AND MRM16 != PM
                    AND MRQ76  != PH9
                    AND MRQ76  != MASRIA
                    AND MRQ76  != MR297
                    AND MRQ76  != MR307
                    AND MRQ76  != MR12H
                    AND MRQ76  != MR127
                    AND MRQ76  != MRQ74
                    AND MRQ76  != MR103
                    AND MRQ76  != MR106
                    AND MRQ76  != PM 
                    AND PS != PH9
                    AND PS != MASRIA
                    AND PS != MR297
                    AND PS != MR307
                    AND PS != MR12H
                    AND PS != MR127
                    AND PS != MRQ74
                    AND PS != MR103
                    AND PS!= MR106
                    AND PS!= PM ") 


medium_pt<- sqldf("select * from hapmap where MR297 == MR307
                    AND MR297 == MR12H
                    AND MR12H == MR307
                    AND MR297 != PH9
                    AND MR297 != MASRIA
                    AND MR297 != MRM16
                    AND MR297 != MRQ76
                    AND MR297 != PS
                    AND MR297 != MR127
                    AND MR297 != MRQ74
                    AND MR297 != MR103
                    AND MR297 != MR106
                    AND MR297 != PM
                    AND MR307 != PH9
                    AND MR307 != MASRIA
                    AND MR307 != MRM16
                    AND MR307 != MRQ76
                    AND MR307 != PS
                    AND MR307 != MR127
                    AND MR307 != MRQ74
                    AND MR307 != MR103
                    AND MR307 != MR106
                    AND MR307 != PM
                    AND MR12H != PH9
                    AND MR12H != MASRIA
                    AND MR12H != MRM16
                    AND MR12H != MRQ76
                    AND MR12H != PS
                    AND MR12H != MR127
                    AND MR12H != MRQ74
                    AND MR12H != MR103
                    AND MR12H != MR106
                    AND MR12H != PM ") 


hard_pt<- sqldf("select * from hapmap where MR127 == MRQ74
                    AND MR127 == MR103
                    AND MR127 == MR106
                    AND MR103 == MR106
                    AND MR103 == MRQ74
                    AND MR106 == MRQ74
                    AND MR127 != 'NC'
                    AND MRQ74 != 'NC'
                    AND MR103 != 'NC'
                    AND MR106 != 'NC'
                    AND MR127 != PH9
                    AND MR127 != MASRIA
                    AND MR127 != MRM16
                    AND MR127 != MRQ76
                    AND MR127 != PS
                    AND MR127 != MR297
                    AND MR127 != MR307
                    AND MR127 != MR12H
                    AND MRQ74 != PH9
                    AND MRQ74 != MASRIA
                    AND MRQ74 != MRM16
                    AND MRQ74  != MRQ76
                    AND MRQ74  != PS
                    AND MRQ74  != MR297
                    AND MRQ74  != MR307
                    AND MRQ74  != MR12H
                    AND MR103 != PH9
                    AND MR103 != MASRIA
                    AND MR103 != MRM16
                    AND MR103  != MRQ76
                    AND MR103  != PS
                    AND MR103  != MR297
                    AND MR103  != MR307
                    AND MR103  != MR12H
                    AND MR106 != PH9
                    AND MR106 != MASRIA
                    AND MR106 != MRM16
                    AND MR106  != MRQ76
                    AND MR106  != PS
                    AND MR106  != MR297
                    AND MR106  != MR307
                    AND MR106  != MR12H ")  

waxy_pt<- select(waxy_pt, snpID, chr, pos, ref, alel,
                PH9, PM, MASRIA, MRM16, MRQ76, PS,
                MR297, MR307, MR12H, MR127, MRQ74, MR103, MR106)

waxy_pt<- merge(waxy_pt, genic_ann, by=c("snpID", "chr", "pos", "ref", "alel"))

pattern_mbkbase<- merge(pattern_13var, mbkbase, by="R498_ID")


#if select from list of 21var

PM<- genotype %>%
  filter(PM != "NC") %>%
  select(snpID, ref, PM)

tsnp<- Reduce(function(x,y) merge(x, y, by = "snpID", all.x = TRUE, all.y = TRUE),
               list(PM, PH9, MASRIA, MRM16, MRQ76, PS, MR297, MR307, MR12H, 
                    MR127, MRQ74, MR103, MR106))

tsnp<- select(tsnp, snpID, PH9, PM, MASRIA, MRM16, MRQ76, PS,
              MR297, MR307, MR12H, MR127, MRQ74, MR103, MR106)

gsnp<- Reduce(function(x,y) merge(x, y, by = c("snpID", "ref", "alel"), all.x = TRUE, all.y = TRUE),
              list(PM, PH9, MASRIA, MRM16, MRQ76, PS, MR297, MR307, MR12H, 
                   MR127, MRQ74, MR103, MR106))

#hapmap format
geno13var$PM<- str_replace(geno13var$PM, "NC", geno13var$ref)
geno13var$PH9<- str_replace(geno13var$PH9, "NC", geno13var$ref)
geno13var$PS<- str_replace(geno13var$PS, "NC", geno13var$ref)
geno13var$MRQ76<- str_replace(geno13var$MRQ76, "NC", geno13var$ref)
geno13var$MRQ74<- str_replace(geno13var$MRQ74, "NC", geno13var$ref)
geno13var$MASRIA<- str_replace(geno13var$MASRIA, "NC", geno13var$ref)
geno13var$MRM16<- str_replace(geno13var$MRM16, "NC", geno13var$ref)
geno13var$MR297<- str_replace(geno13var$MR297, "NC", geno13var$ref)
geno13var$MR307<- str_replace(geno13var$MR307, "NC", geno13var$ref)
geno13var$MR12H<- str_replace(geno13var$MR12H, "NC", geno13var$ref)
geno13var$MR127<- str_replace(geno13var$MR127, "NC", geno13var$ref)
geno13var$MR103<- str_replace(geno13var$MR103, "NC", geno13var$ref)
geno13var$MR106<- str_replace(geno13var$MR106, "NC", geno13var$ref)

#order chr
a<- filter(hapmap13var, chr %in% c("Chr1", "Chr2", "Chr3", "Chr4", "Chr5", "Chr6", 
                              "Chr7", "Chr8",
                              "Chr9"))
b<- filter(hapmap13var, chr %in% c("Chr10", "Chr11", "Chr12"))

a<- a[order(a$chr, a$pos),]
b<- b[order(b$chr, b$pos),]

k<- rbind(a,b)

#qtl analysis

#6900000-8200000 qSAC3
#Os02g0513000 qAC2   OsR498G0203912300.01
#Os08g0122000 qAC8   OsR498G0815094200.01

qSAC3_pattern<- a %>%
  #filter(chr == "Chr8") %>%
  filter(dplyr::between(pos, 6900000, 8200000)) 

py_qSAC3<- ann_qSAC3 %>%
  filter(str_detect(alel, "/"))

#extract based on string
shatter <- shatter %>% 
  filter(str_detect(R498_ID, "OsR498G"))

sum_qSAC3<- sqldf("select count(snpID) as countSNP,  R498_ID, geneDesc
               from ann_qSAC3 group by R498_ID")

a<- a %>% 
  mutate(across(everything(), as.character))#all as.character


#nest & paste
k<- kegg_ceri %>%
  group_by(geneID) %>%
  summarise(Pathway.ID=paste(Pathway.ID, collapse = ";")) 
