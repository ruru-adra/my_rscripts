#sequenom conversion
s3<- s3 %>% mutate(across(everything(), as.character)) #change to char

s3<- unique(s3) # remove dup

s3$alel[s3$alel==''] <- NA #replace missing value with NA

s3$alel[is.na(s3$alel)]<- "NC" #replace NA with string value

n2<- n %>% #convert long to wider
  pivot_wider(names_from = sampleID, values_from = alel) %>% 
  unnest(cols = everything() )

*********

library(tidyr)
library(dplyr)
library(stringi)
library(stringr)
library(sqldf)
library(CMplot)
library(readxl)


raw_PH9<- read.table("recode_INDEL_PH9_gatk.txt", header=T, sep="\t", quote="")
raw_MSRIA<- read.table("INDEL_cc_MASRIA.recode.txt", header=T, sep="\t", quote="")
raw_MRM16<- read.table("recode_INDEL_MRM16_gatk.txt", header=T, sep="\t", quote="")
raw_Q76<- read.table("recode_INDEL_Q76_gatk.txt", header=T, sep="\t", quote="")
raw_PS<- read.table("recode_INDEL_PS_gatk.txt", header=T, sep="\t", quote="")
raw_MR297<- read.table("recode_INDEL_MR297_gatk.txt", header=T, sep="\t", quote="")
raw_MR307<- read.table("INDEL_cc_MR307.recode.txt", header=T, sep="\t", quote="")
raw_MR12H<- read.table("../../../../HYBRID/fltr_var/recode_cc_KADARIA_indel.txt", header=T, sep="\t", quote="")
raw_MR127<- read.table("recode_INDEL_MR127_gatk.txt", header=T, sep="\t", quote="")
raw_Q74<- read.table("INDEL_cc_Q74.vcf.recode.txt", header=T, sep="\t", quote="")
raw_MR103<- read.table("INDEL_cc_MR103.recode.txt", header=T, sep="\t", quote="")
raw_MR106<- read.table("INDEL_cc_MR106.recode.txt", header=T, sep="\t", quote="")
raw_PM<- read.table("Results/VARIANT/SNP_PM_fltr_all.txt", header=T, sep="\t", quote="")

all<- read.table("Results/VARIANT/genotype_snp_cohort_amylose_latest.txt", header=T, sep="\t", quote="")

#colnames
names(raw_PM)<- c("chr", "ref", "alel", "qual", "pos", "GT", "AD", "DP", "GQ", "PL")

#filter alel & DP
fltr_PM<- raw_PM %>% filter(DP >= 2)

fltr_PM$snpID<- paste(fltr_PM$chr, fltr_PM$pos, sep="_")

idl_MR103$indelID<- paste(idl_MR103$chr, idl_MR103$pos, sep="_")
idl_MR106$indelID<- paste(idl_MR106$chr, idl_MR106$pos, sep="_")
idl_MR127$indelID<- paste(idl_MR127$chr, idl_MR127$pos, sep="_")
idl_MR12H$indelID<- paste(idl_MR12H$chr, idl_MR12H$pos, sep="_")
idl_MR297$indelID<- paste(idl_MR297$chr, idl_MR297$pos, sep="_")
idl_MR307$indelID<- paste(idl_MR307$chr, idl_MR307$pos, sep="_")
idl_MRM16$indelID<- paste(idl_MRM16$chr, idl_MRM16$pos, sep="_")
idl_MSRIA$indelID<- paste(idl_MSRIA$chr, idl_MSRIA$pos, sep="_")
idl_PH9$indelID<- paste(idl_PH9$chr, idl_PH9$pos, sep="_")
idl_PS$indelID<- paste(idl_PS$chr, idl_PS$pos, sep="_")

#snp in genotype format
g_Q76<- unique(select(Q76, snpID, ref, Q76))
 
y_snp<- Reduce(function(x,y) merge(x, y, by = c("snpID", "ref", "alel"), all.x = TRUE, all.y = TRUE),
             list(g_PM, g_PH9, g_MASRIA, g_MRM16, g_Q76, g_PS, g_MR297, g_MR307, g_MR12H, g_MR127,
                  g_Q74, g_MR103, g_MR106))

g_snp<- Reduce(function(x,y) merge(x, y, by = "snpID", all.x = TRUE, all.y = TRUE),
                 list(g_PM, g_PH9, g_MASRIA, g_MRM16, g_Q76, g_PS, g_MR297, g_MR307, g_MR12H, g_MR127,
                      g_Q74, g_MR103, g_MR106))


snp<- select(g_snp, snpID, PM, PH9, MASRIA, MRM16, Q76, PS, MR297, MR307, MR12H, MR127, Q74, MR103, MR106)

all2<- unique(select(all, snpID, ref, alel))

a<- sqldf("select count(indelID) as count_indel, indelID from t_indel group by indelID")
a<- sqldf("select count(snpID) as count_snp, snpID from y_snp group by snpID")

indel_var<- filter(a, count_indel != 1)
indel_sgle<- filter(a, count_indel == 1)
snp_var<- filter(a, count_snp != 1)

indel_sgle<- indel_sgle %>%
  merge(t_indel, by="indelID") %>%
  merge(indel, by="indelID")

indel_var2<- indel_var %>%
  merge(t_indel, by="indelID") %>%
  merge(indel, by="indelID")

snp_var2<- snp_var %>%
  merge(y_snp, by="snpID") %>%
  merge(indel, by="indelID")

ref_sama<- sqldf("select indelID, count_indel, ref, alel from indel_var2 t1
                 WHERE EXISTS (select indelID, ref, alel from indel_var2 t2
                 WHERE t1.indelID == t2.indelID AND
                       t1.ref == t2.ref)") #if sama/taksama replace ==/<>

ref_sama<- anti_join(indel_var2, ref_taksama, by="indelID") #polymorphic


ref<- tsnp[, c(1,2,4,6,8,10,12,14,16,18,20,22,24)]
names(ref)<- c("snpID", "ref_PH9", "ref_MASRIA", "ref_MRM16", "ref_Q76", "ref_PS", "ref_MR297", 
          "ref_MR307", "ref_MR12H", "ref_MR127", "ref_Q74", "ref_MR103", "ref_MR106")

a12<- ref %>%
  select(snpID, ref_MR106) %>%
  filter(ref_MR106 != "NA") %>%
  unique()
colnames(a12)[2]="ref"


ref_snp<- unique(rbind(a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12))


my_snp<- merge(ref_snp, snp, by="snpID")

tsnp$PH9<- as.character(tsnp$PH9)
tsnp$MASRIA<-  as.character(tsnp$MASRIA)
tsnp$MRM16<-  as.character(tsnp$MRM16) 
tsnp$Q76<-  as.character(tsnp$Q76) 
tsnp$PS<-  as.character(tsnp$PS) 
tsnp$MR297<-  as.character(tsnp$MR297) 
tsnp$MR307<-  as.character(tsnp$MR307) 
tsnp$MR12H<-  as.character(tsnp$MR12H) 
tsnp$MR127<-  as.character(tsnp$MR127)
tsnp$Q74<-  as.character(tsnp$Q74) 
tsnp$MR103<-  as.character(tsnp$MR103) 
tsnp$MR106<-  as.character(tsnp$MR106)

tsnp[is.na(tsnp)]<- "NC"

tsnp<- tsnp[, c(1,3,4,5,6,7,8,9,10,11,12,13,14)]

tsnp$count<- rowSums(tsnp == "NC")
tsnp$A<- rowSums(tsnp == "A")
tsnp$C<- rowSums(tsnp == "C")
tsnp$G<- rowSums(tsnp == "G")
tsnp$T<- rowSums(tsnp == "T")

pv_alel<- filter(tsnp, count == "11") #unique snp

py_count<- select(tsnp, snpID, A, C, G, T)
py_count$countZero<- rowSums(py_count == 0)
py_alel<- filter(py_count, countZero != 3)

mono_alel<- anti_join(tsnp, a_b, by="snpID")

#result py
#Chr1 Chr10 Chr11 Chr12  Chr2  Chr3  Chr4  Chr5  Chr6  Chr7  Chr8  Chr9 
#1078   444  1129   668  1163   473  1290   285   411   278   882   576
#rate py - 0.28%

tsnp_my<- merge(my_snp, chr_snp, by="snpID")

m<- tsnp %>%
  filter(count == 0) %>%
  select(snpID, count)

write.table(fltr_snp, "list_fltr_snp_amylose.txt", col.names=T, row.names=F, sep="\t", quote=F)
write.table(tsnp_my, "snp_ngs_amylose.txt", col.names=T, row.names=F, sep="\t", quote=F)
write.table(my_genotype, "genotype_my.txt", col.names=T, row.names=F, sep="\t", quote=F)

#create alel
colA<- py_count %>%
  filter(A != 0) %>%
  select(snpID)

colC<- py_count %>%
  filter(C != 0) %>%
  select(snpID)

colG<- py_count %>%
  filter(G != 0) %>%
  select(snpID)

colT<- py_count %>%
  filter(T != 0) %>%
  select(snpID)

colA$alel="A"
colC$alel="C"
colG$alel="G"
colT$alel="T"

alel<- rbind(colA, colC, colG, colT)
alel_py<- alel %>%
  merge (tsnp_my, by="snpID") %>%
  filter(type=="polymorphic")
  
#create alel for polymorphic snp
new_a<- a %>%
     pivot_wider(names_from = countZero, values_from= alel)

my_genotype<- merge(tsnp_my, alel, by="snpID")

#for indel create refID #as.character

new_a<- a %>%
  pivot_wider(names_from = ref, values_from= alel)

aa<- a2 %>%
  filter(A != "") %>%
  select(indelID, A) %>%
  unique()

py_indel<- rbind(aa, ac, ag, at)

indel_multiref<- ref_taksama %>% 
                 merge(indel_var, by=c("indelID", "ref", "alel")) %>%
                 select(indelID, ref, alel, PH9, MASRIA, MRM16, Q76, PS, MR297, MR307, MR12H,
                        MR127, Q74, MR103, MR106) %>%
                 unique()

snp <- snp[order(snp$chr, snp$pos),] #pos as.numeric

#annotation

vep<- read.table("vep_to_vcf_cohort_indel.vcf", header=T, sep="\t", quote="", skip=28 )
seff<- read.table("ann_indel_amylose_ngs.txt", header=T, sep="\t", quote="" )

names(seff)<- c("chr", "pos", "ref", "alel", "effect", "OsID", "AA", "AA_POS", "IMPACT", "RANK")

fltr_seff<- seff %>%
  filter(effect %in% c("intergenic_region", "missense_variant", "intron_variant",
                       "synonymous_variant", "5_prime_UTR_variant", "3_prime_UTR_variant")) %>%
  select(chr, pos, ref, alel, effect, OsID, AA, AA_POS, IMPACT, RANK) %>%
  unique()

fltr_seff$snpID<- paste(fltr_seff$chr, fltr_seff$pos, sep="_")
fltr_seff$indelID<- paste(fltr_seff$chr, fltr_seff$pos, sep="_")

fltr_seff2<- unique(select(fltr_seff, indelID, chr, pos, ref, alel, effect, OsID))
  
gby_seff<- sqldf("select count(effect) as count_effect, indelID, effect from fltr_seff2
                 group by indelID")
gby_seff<- arrange(gby_seff, desc(count_effect))

sgl_ann<- filter(gby_seff, count_effect == 1)
dup_ann<- filter(gby_seff, count_effect != 1)

dup_seff<- merge(dup_ann, fltr_seff2, by="indelID")

sgl_seff<- merge(sgl_ann, fltr_seff2, by="indelID")

#vep data clean
vep<- read.table("../ann_var/vep_amylose.vcf", header=F, sep="\t", quote="", skip=28)

vep<- vep[, c(1,4,7,10,11,14)]

names(vep)<- c("indelID", "OsID", "effect", "AA_POS", "AA", "IMPACT")

vep<- separate(data=vep, col=indelID, into=c("chr", "pos", "ref_alel"), sep="_")
vep<- separate(data=vep, col=ref_alel, into=c("ref", "alel"), sep="/")
vep<- separate(data=vep, col=IMPACT, into="IMPACT", sep=";")

fltr_vep<- vep %>%
  filter(effect %in% c("intergenic_region", "missense_variant", "intron_variant",
                       "synonymous_variant", "5_prime_UTR_variant", "3_prime_UTR_variant")) %>%
  select(chr, pos, ref, alel, effect, OsID, AA, AA_POS, IMPACT) %>%
  unique()

fltr_vep$snpID<- paste(fltr_vep$chr, fltr_vep$pos, sep="_")
fltr_vep$indelID<- paste(fltr_vep$chr, fltr_vep$pos, sep="_")

fltr_vep2<- unique(select(fltr_vep, snpID, chr, pos, ref, alel, effect, OsID))
fltr_vep2<- unique(select(fltr_vep, indelID, chr, pos, ref, alel, effect, OsID))

#visualisation

to_plot<- select(snp, snpID, chr, pos)

CMplot(ab,type="p",plot.type="d",bin.size=1e6,chr.den.col=c("darkgreen", "yellow", "red"),file="jpg",memo="",dpi=300,
       file.output=TRUE,verbose=TRUE,width=9,height=6)

#mkbase
shatter<- read_excel("~/Downloads/mkbase_quality.xlsx", sheet=11)
names(shatter)<- c("geneSymbol", "LocID", "R498_ID", "Trait", "Class")
shatter<- shatter %>%
  mutate(R498_ID=strsplit(as.character(R498_ID),"\\s+")) %>%
  unnest(R498_ID)
shatter$R498_ID<- str_trim(shatter$R498_ID)
shatter <- shatter %>% 
  filter(str_detect(R498_ID, "OsR498G"))
shatter$Trait="seed shattering"
shatter$Class="Quality"
shatter<- R498 %>%
  merge(shatter, by="R498_ID") %>%
  select(R498_ID, geneSymbol, geneDesc, Trait, Class) %>%
  unique()

#annotation

ag_snp<- ann_snp %>%
  merge(ag_trait, by="R498_ID") %>%
  select(snpID, chr, pos, ref, alel, R498_ID, type, effect, geneSymbol, geneDesc, Trait, Class) %>%
  unique()

ag_gby<- sqldf("select count(distinct R498_ID) as count_gene, snpID
               from ag_snp group by snpID") #count how many gene per snp 
                                            #total snp with ag annotation=35,313
ag_gby<- arrange(ag_gby, desc(count_gene))


snp_per_gene<- sqldf("select count(distinct snpID) as count_snp, R498_ID
               from ag_snp group by R498_ID") # 1669 genes containing snps
snp_per_gene<- arrange(snp_per_gene, desc(count_snp))

snp_per_gene2<- merge(snp_per_gene, ag_trait, by="R498_ID")


snp_quality<- filter(snp_per_gene2, Class=="Quality")

uq<- unique(select(snp_morpho, R498_ID))
ts<- unique(select(snp_morpho, R498_ID, count_snp))

write.table(snp_abiotic, "./Results/abiotic_snp.txt", col.names=T, row.names=F, sep="\t", quote=F)

py_snp<- filter(ag_snp, type=="polymorphic")

nsyn_per_gene<- sqldf("select count(distinct snpID) as count_snp, R498_ID
               from ns_snp group by R498_ID")

#amylose snp

snp_ssrg<- merge(ann_snp, ssrg, by="R498_ID")

nsyn_ssrg<- filter(snp_ssrg, effect=="missense_variant")

gby_effect<- sqldf("select count(distinct snpID) as count_snp, effect
               from snp_ssrg group by effect")

gby_cat<- sqldf("select count(distinct snpID) as count_snp, category
               from snp_ssrg group by category")

gby_var<- sqldf("select count(distinct snpID) as count_snp, category
               from snp_ssrg group by ")


cd2tassel<- cd %>%
  filter(GROUP %in% c("H", "S", "W", "M")) %>%
  select(snpID) %>%
  unique()


hard<- snp_ssrg %>%
  filter(MR127 != "NA") %>%
  filter(Q74 != "NA") %>%
  filter(MR103 != "NA") %>%
  filter(MR106 != "NA") %>%
  select(R498_ID, snpID, chr, pos, ref, alel, MR127, Q74, MR103, MR106) %>%
  unique()


##########################################################

MR106<- all %>%
  filter(MR106 != "NC") %>%
  select(snpID, chr, pos, ref, alel, MR106) %>%
  unique()


snp<- select(g_snp, snpID, PM, PH9, MASRIA, MRM16, Q76, PS, MR297, MR307, MR12H, MR127, Q74, MR103, MR106)

all2<- unique(select(all, snpID, ref, alel))

snp_p2<- merge(snp, all2, by="snpID")

mb<- merge(m, b, by="snpID")

a$PM<- as.character(snp$PM)
snp$PH9<- as.character(snp$PH9)
snp$MASRIA<-  as.character(snp$MASRIA)
snp$MRM16<-  as.character(snp$MRM16) 
snp$Q76<-  as.character(snp$Q76) 
snp$PS<-  as.character(snp$PS) 
snp$MR297<-  as.character(snp$MR297) 
snp$MR307<-  as.character(snp$MR307) 
snp$MR12H<-  as.character(snp$MR12H) 
snp$MR127<-  as.character(snp$MR127)
snp$Q74<-  as.character(snp$Q74) 
snp$MR103<-  as.character(snp$MR103) 
snp$MR106<-  as.character(snp$MR106)

snp[is.na(snp)]<- "NC"

snp<- filter(snp, chr %in% c("Chr1", "Chr2", "Chr3", "Chr4", "Chr5", "Chr6", "Chr7", "Chr8",
                             "Chr9", "Chr10", "Chr11", "Chr12"))
                             
                             x<- filter(cnd2primer, chr %in% c("Chr1", "Chr2", "Chr3", "Chr4", "Chr5", "Chr6", "Chr7", "Chr8",
                             "Chr9"))

tsnp<- tsnp[, c(1,3,4,5,6,7,8,9,10,11,12,13,14)]

tsnp$countNC<- rowSums(tsnp == "NC")
tsnp$A<- rowSums(tsnp == "A")
tsnp$C<- rowSums(tsnp == "C")
tsnp$G<- rowSums(tsnp == "G")
tsnp$T<- rowSums(tsnp == "T")

pv_alel<- filter(tsnp, countNC == "12") #unique snp

py_count<- select(tsnp, snpID, A, C, G, T)
py_count$countZero<- rowSums(py_count == 0)
py_alel<- filter(py_count, countZero != 3)

mono_alel<- anti_join(tsnp, ab, by="snpID")



amy_snp<- merge(snp, char_snp, by="snpID")

m<- tsnp %>%
  filter(count == 0) %>%
  select(snpID, count)

new_a<- snp_var2 %>%
  pivot_wider(names_from = head, values_from= alel)


###amylose ann

names(ann)<- c("chr", "pos", "ref", "alel", "effect", "OsID", "AA", "AA_POS", "IMPACT", "RANK")

fltr_ann<- ann %>%
  filter(effect %in% c("intergenic_region", "missense_variant", "intron_variant",
                       "synonymous_variant", "5_prime_UTR_variant", "3_prime_UTR_variant")) %>%
  select(chr, pos, ref, alel, effect, R498_ID, AA, AA_POS, IMPACT, RANK) %>%
  unique()

fltr_ann$snpID<- paste(fltr_ann$chr, fltr_ann$pos, sep="_")#3,113,000 annotated snps


no_ann<- anti_join(amy_snp, uq, by="snpID")

nsyn_13var<- filter(fltr_ann, effect=="missense_variant")
















