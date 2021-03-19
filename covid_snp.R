
cov2_snp<- read.table("../MGI/COV2-1_SNP.txt", header=T, sep="\t", quote="")
cov2_indel<- read.table("../MGI/COV2-1.INDEL.txt", header=T, sep="\t", quote="")

cov1_snp$snpID<- paste("snp", cov1_snp$POS, sep="_")#create snpID
cov2_snp$snpID<- paste("snp", cov2_snp$POS, sep="_")
cov3_snp$snpID<- paste("snp", cov3_snp$POS, sep="_")
cov4_snp$snpID<- paste("snp", cov4_snp$POS, sep="_")
cov5_snp$snpID<- paste("snp", cov5_snp$POS, sep="_")
cov6_snp$snpID<- paste("snp", cov6_snp$POS, sep="_")

s_cov1<- unique(select(cov1_snp, snpID, REF, ALT)) #select alel for snp table
s_cov2<- unique(select(cov2_snp, snpID, REF, ALT))
s_cov3<- unique(select(cov3_snp, snpID, REF, ALT))
s_cov4<- unique(select(cov4_snp, snpID, REF, ALT))
s_cov5<- unique(select(cov5_snp, snpID, REF, ALT))
s_cov6<- unique(select(cov6_snp, snpID, REF, ALT))

tsnp<- Reduce(function(x,y) merge(x, y, by = c("snpID", "REF"), all.x = TRUE, all.y = TRUE), #run snp table
              list(s_cov1, s_cov2, s_cov3, s_cov4, s_cov5, s_cov6))

tsnp$COV1<- as.character(tsnp$COV1)
tsnp$COV2<- as.character(tsnp$COV2)
tsnp$COV3<- as.character(tsnp$COV3)
tsnp$COV4<- as.character(tsnp$COV4)
tsnp$COV5<- as.character(tsnp$COV5)
tsnp$COV6<- as.character(tsnp$COV6)

tsnp[is.na(tsnp)]<- "NC"


#filter position in genic region
snp_orf1ab<- tsnp %>%
  filter(dplyr::between(pos, 266, 21555))

snp_spike<- tsnp %>%
  filter(dplyr::between(pos, 21563, 25384))

snp_ORF3a<- tsnp %>%
  filter(dplyr::between(pos, 25393,	26220))

snp_envelope<- tsnp %>%
  filter(dplyr::between(pos, 26245, 26472))

snp_mglyco<- tsnp %>%
  filter(dplyr::between(pos, 26523,	27191))

snp_ORF6<- tsnp %>%
  filter(dplyr::between(pos, 26523,	27191))

snp_ORF7a<- tsnp %>%
  filter(dplyr::between(pos, 26523,	27191))

snp_ORF6<- tsnp %>%
  filter(dplyr::between(pos, 27202,	27387))

snp_ORF7a<- tsnp %>%
  filter(dplyr::between(pos, 27394,	27759))

snp_ORF7b<- tsnp %>%
  filter(dplyr::between(pos, 27756,	27887))

snp_ORF8<- tsnp %>%
  filter(dplyr::between(pos, 27894,	28259))

snp_nucleo<- tsnp %>%
  filter(dplyr::between(pos, 28274,	29533))

snp_ORF10<- tsnp %>%
  filter(dplyr::between(pos, 29558,	29674))



#266	13468	CDS	product=orf1ab polyprotein	YP_009724389.1
#13468	21555	CDS	product=orf1ab polyprotein	YP_009724389.1
#266	13483	CDS	product=orf1a polyprotein	YP_009725295.1
#21563	25384	CDS	product=surface glycoprotein	YP_009724390.1
#25393	26220	CDS	product=ORF3a protein	YP_009724391.1
#26245	26472	CDS	product=envelope protein	YP_009724392.1
#26523	27191	CDS	product=membrane glycoprotein	YP_009724393.1
#27202	27387	CDS	product=ORF6 protein	YP_009724394.1
#27394	27759	CDS	product=ORF7a protein	YP_009724395.1
#27756	27887	CDS	product=ORF7b	YP_009725296.1
#27894	28259	CDS	product=ORF8 protein	YP_009724396.1
#28274	29533	CDS	product=nucleocapsid phosphoprotein	YP_009724397.2
#29558	29674	CDS	product=ORF10 protein	YP_009725255.1
#1	265	five_prime_UTR		
#29675	29903	three_prime_UTR		

snp_5UTR<- tsnp %>%
  filter(dplyr::between(pos, 1, 265))

snp_3UTR<- tsnp %>%
  filter(dplyr::between(pos, 29675,	29903))

snp_14408<- filter(msa, snpID=="snp_14408")
snp_snp_17747<- filter(msa, snpID=="snp_17747")
snp_17858<- filter(msa, snpID=="snp_17858")
snp_18060<- filter(msa, snpID=="snp_18060")
snp_23403<- filter(msa, snpID=="snp_23403")
snp_28144<- filter(msa, snpID=="snp_28144")
snp_3037<- filter(msa, snpID=="snp_3037")
snp_8782<- filter(msa, snpID=="snp_8782")
snp_241<- filter(msa, snpID=="snp_241")

wide_14408<-  snp_14408 %>% 
  spread(sample, gtype)

wide_241<-  snp_241 %>% 
  spread(sample, gtype)

wide_17858<-  snp_17858 %>% 
  spread(sample, gtype)

wide_18060<-  snp_18060 %>% 
  spread(sample, gtype)

wide_23403<-  snp_23403 %>% 
  spread(sample, gtype)

wide_28144<-  snp_28144 %>% 
  spread(sample, gtype)

wide_3037 <-  snp_3037 %>% 
  spread(sample, gtype)

wide_8782<-  snp_8782 %>% 
  spread(sample, gtype)

wide_17747<-  snp_17747 %>% 
  spread(sample, gtype)

a<- rbind(wide_14408, wide_17747, wide_17858, wide_18060, wide_23403,
          wide_241, wide_3037, wide_8782, wide_28144)

wide_snp_msia<-  snp_msia %>% 
  spread(sample_name, allele)


wide<-  seq1009 %>% 
  spread(sampleID, ref_allele)