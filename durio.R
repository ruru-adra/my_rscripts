#read
dz_snp<- read.table("dz_durio_snp.txt", header=T, sep="\t", quote="")
names(dz_snp)<- c("chr", "ref", "alel", "qual", "pos", "GT", "AD",
                 "DP", "GQ", "PL")


fltr_l<- l_snp %>%
  filter(alel %in% c("A", "C", "G", "T"))

fltr_l$snpID<- paste(fltr_l$chr, fltr_l$pos, sep="_")
  

g_dz<- unique(select(fltr_dz, snpID, ref, alel))
g_g<- unique(select(fltr_g, snpID, ref, alel))
g_k<- unique(select(fltr_k, snpID, ref, alel))
g_l<- unique(select(fltr_l, snpID, ref, alel))


colnames(g_dz)[3]="alel"
colnames(g_g)[3]="alel"
colnames(g_k)[3]="alel"
colnames(g_l)[3]="alel"

tsnp<- Reduce(function(x,y) merge(x, y, by = "snpID", all.x = TRUE, all.y = TRUE),
              list(g_dz, g_g, g_k, g_l))

gsnp<- Reduce(function(x,y) merge(x, y, by = c("snpID", "ref", "alel"), all.x = TRUE, all.y = TRUE),
              list(g_dz, g_g, g_k, g_l)) #get py allele

n <- n %>%
  mutate_if(sapply(n, is.factor), as.character) # all col

a<- tsnp
a[is.na(a)]<- "NC"
a$countNC<- rowSums(a == "NC")

uq_snp<- filter(a, countNC== 3)#unique

t<- filter(a, countNC==0)

mono_snp<- t #mono
prob_py<- merge(mono_snp, d_snp, by="snpID")

py<- anti_join(d_snp, prob_py, by="snpID")
py<- merge(py, tsnp, by="snpID")

n$DZ<- str_replace(n$DZ, "NC", n$ref)
n$G<- str_replace(n$L, "NC", n$ref)
n$K<- str_replace(n$K, "NC", n$ref)
n$L<- str_replace(n$L, "NC", n$ref)

#merge py + prob py

n2<- gsnp %>% #pivot to get py allele in single row
  pivot_wider(names_from = alel, values_from=alel)
n2$alel<- paste(n2$G, n2$A, n2$T, n2$C, sep="/")




names(d)<- c("chr", "pos", "dot", "ref", "alel", "dot1", "dot2", "dot3")
d<- d[, c(1,1,2,3,4,5,6,7,8)]
colnames(d)[2]="chr_new"
d<- separate(data=d, col=chr_new, into=c("n", "id"), sep="_")

d$pos<- as.numeric(d$pos)
d$id<- as.numeric(d$id)
d2<- d[order(d$id, d$pos),]



utr<- genic %>%
  filter(effect %in% c("3_prime_UTR_variant","5_prime_UTR_variant")) %>%
  select(snpID) %>%
  unique()

dt$pos_diff <- c(0, diff(dt$pos))












names(ann_a)<- c("chr", "pos", "ref", "alel", "eff", "R498_ID", "AA", "POS", "IMPACT")

a<- filter(ann, eff %in% c("intergenic_region", "missense_variant", 
                           "synonymous_variant", "intron_variant"))

h_tsl$strand="+"
h_tsl$assembly="."
h_tsl$assembly="HYBRID"
h_tsl$protLSID="NA"
h_tsl$assayLSID="NA"
h_tsl$panel="hybrid"
h_tsl$QCcode="NA"
h_tsl$center="BN"



NBS_IPR002182<- fltr_ips %>%
  filter(str_detect(IProID, "IPR002182"))

TIR_IPR000157<- fltr_ips %>%
  filter(str_detect(IProID, "IPR000157"))

LRR1_IPR001611<- fltr_ips %>%
  filter(str_detect(IProID, "IPR001611"))

LRR2_IPR013101<- fltr_ips %>%
  filter(str_detect(IProID, "IPR013101"))

LRR3_IPR011713<- fltr_ips %>%
  filter(str_detect(IProID, "IPR011713"))

to_pcc<- fpkm %>% pivot_wider(names_from= GeneID, values_from = log2FoldChange)


clean<- fpkm %>%
filter(BaliBR != "0") %>%
  filter(MR297WR != "0") %>%
  filter(MRM16RR != "0") %>%
  filter(MRQ100RR != "0") %>%
  filter(MRQ76WR != "0") %>%
  filter(PH9BR != "0")

a$cor<- round(a$cor, digits=2)



a %>% filter(between(corr, 0.80, 0.89) )

rmk11_matdays<- rmk11 %>%
  merge(matdays, by="acc_num") %>%
  filter(mat_days <= 150)

amy_mat <- amylose %>%
  merge(matdays, by="acc_num") %>%
  filter(between(mat_days, 0.1, 150)) %>%
  filter(Amylose != 0)

a1<- pcc_fltr %>% filter(between(cor,0.01,0.09))

clstr<- clstr %>%
  mutate(OsID= strsplit(as.character(OsID), ",")) %>%
  unnest(OsID)

japo_flank$bkt<- "["
japo_flank$bkt2<- "]"
japo_flank$bkt3<- "/"

japo_flank$seq<- paste(japo_flank$bkt, japo_flank$ref, japo_flank$bkt3, japo_flank$alel, japo_flank$bkt2, sep="")

japo_flank$flank<- paste(japo_flank$left_seq, japo_flank$seq, japo_flank$right_seq, sep="")

a2<- a %>%
  merge(cd, by="snpID") %>%
  select(snpID, chr.x, ref.x, alel.x, flank, PH9, PM, MASRIA, MRM16,  MRQ76,    
          PS, MR297,MR307,MR12H,MR127,MRQ74,MR103,MR106)

b2<- b %>%
  merge(chr6, by="snpID") %>%
  select(snpID, chr.x, ref.x, alel.x, flank, PH9, PS, MRQ76, MR127)

a<- mcl %>%
  mutate(OsID= strsplit(as.character(OsID), ";")) %>%
  unnest(OsID)

mcl_C10<- mcl %>%
  filter(clstrID=="C10") %>%
  merge(rapdb,by="OsID")

mcl_go<- mcl_go %>%
  mutate(CC= strsplit(as.character(CC), ";")) %>%
  unnest(CC)

var_ssrg<- filter(var_ssrg, var %in% c("MASRIA", "MR103", "MR106", "MR127", "MR12H", "MR297", "MR307",
                                 "MRM16", "PH9", "PM", "PS", "MRQ74", "MRQ76"))

a<- filter(a, ALT %in% c("A", "C", "G", "T"))

a<- geno_up %>%
  pivot_wider(names_from = var, values_from=ALT)


a<- a %>% 
  mutate(across(everything(), as.character))#all as.character

a[is.na(a)]<- "NC"



a$PM<- str_replace(a$PM, "NC", a$REF)

a$PH9<- str_replace(a$PH9, "NC", a$REF)
a$PS<- str_replace(a$PS, "NC", a$REF)
a$MRQ76<- str_replace(a$MRQ76, "NC", a$REF)
a$MRQ74<- str_replace(a$MRQ74, "NC", a$REF)
a$MASRIA<- str_replace(a$MASRIA, "NC", a$REF)
a$MRM16<- str_replace(a$MRM16, "NC", a$REF)
a$MR297<- str_replace(a$MR297, "NC", a$REF)
a$MR307<- str_replace(a$MR307, "NC", a$REF)
a$MR12H<- str_replace(a$MR12H, "NC", a$REF)
a$MR127<- str_replace(a$MR127, "NC", a$REF)
a$MR103<- str_replace(a$MR103, "NC", a$REF)
a$MR106<- str_replace(a$MR106, "NC", a$REF)


kegg2<- kegg %>%
  mutate(protID= strsplit(as.character(protID), ";")) %>%
  unnest(protID)


ggplot(a, aes(x=GO.Name, y=P.Value)) + 
  geom_bar(stat = "identity", fill="#69b3a2", color="#e9ecef") +
  coord_flip()

ggplot(ab) +
  geom_point(aes(x = P.Value, y= GO.Name, color = GO.Category, alpha=0.5))
               

fltr_react %>%
  arrange(count_seq) %>%    # First sort by val. This sort the dataframe but NOT the factor levels
  #mutate(pathwayName=factor(pathwayName, levels=pathwayName)) %>%   # This trick update the factor levels
  ggplot( aes(x=pathwayName, y=count_seq)) +
  geom_segment( aes(xend=pathwayName, yend=0)) +
  geom_point( size=4, color="green") +
  coord_flip() +
  #theme_bw()
  theme(axis.text.y = element_text(size=12, 
                                    color="black", 
                                    face="bold",
                                    angle=0))
  xlab("")


a<- fltr_react %>%
  mutate(seq= strsplit(as.character(seq), ";")) %>%
  unnest(seq)

