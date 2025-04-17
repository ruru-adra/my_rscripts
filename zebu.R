./snphylo.sh -H /home/biotech/Desktop/snp_indel/snphylo/zebu_5breeds_to_tree.txt -A -t 10

# Select the genotype columns (assuming columns 4 onward)
geno_cols <- names(dt)[4:ncol(dt)]

dt[, simplified_allele := apply(.SD, 1, function(x) {
  paste0(sort(unique(x)), collapse = "/")
}), .SDcols = geno_cols]

setorder(m, chr, pos)

#from wide geno by breed to long geno
a<- hapmap_5breeds %>% 
pivot_longer(cols = BRAHMAN:KKE,
               names_to = "var", 
               values_to = "allele")

#summarize allele per snp
a2<- a %>%
  group_by(snpID) %>%
  summarise(allele=paste(allele, collapse = "/"))

#summarize allele per snp using data.table
dt[, simplified_allele := {
alleles<- tstrsplit(allele, "/", fixed = TRUE)
bases<- unique(unlist(strsplit(unlist(alleles), "")))
paste0(sort(bases), collapse = "/")
}, by = 1:nrow(dt)]

library(dplyr)
library(tidyr)
library(tidyverse)
library(stringi)
library(stringr)


#clean snp
colnames(sh)<- c("chr", "ref", "alel", "qual", "pos", "GT", "AD", "DP", "GQ", "PL")
colnames(lq)<- c("chr", "ref", "alel", "qual", "pos", "GT", "AD", "DP", "GQ", "PL")
colnames(mir)<- c("chr", "ref", "alel", "qual", "pos", "GT", "AD", "DP", "GQ", "PL")
colnames(rcc)<- c("chr", "ref", "alel", "qual", "pos", "GT", "AD", "DP", "GQ", "PL")

sh$snpID <- paste(sh$chr, sh$pos, sep = "_")
lq$snpID <- paste(lq$chr, lq$pos, sep = "_")
mir$snpID <- paste(mir$chr, mir$pos, sep = "_")
rcc$snpID <- paste(rcc$chr, rcc$pos, sep = "_")

fltr_sh<- sh %>%
  filter(alel %in% c("A", "C", "G", "T")) %>%
  filter(DP >= 3) %>%
  filter(qual >= 40) %>%
  filter(GQ >= 30)

fltr_lq<- lq %>%
 filter(alel %in% c("A", "C", "G", "T")) %>%
 filter(DP >= 3) %>%
 filter(qual >= 40) %>%
 filter(GQ >= 30)

fltr_mir<- mir %>%
  filter(alel %in% c("A", "C", "G", "T")) %>%
  filter(DP >= 3) %>%
  filter(qual >= 40) %>%
  filter(GQ >= 30)

fltr_rcc<- rcc %>%
  filter(alel %in% c("A", "C", "G", "T")) %>%
  filter(DP >= 3) %>%
  filter(qual >= 40) %>%
  filter(GQ >= 30)

write.table(fltr_lq, "fltr_lq.txt", col.names=T, row.names=F, sep="\t", quote=F)
write.table(fltr_mir, "fltr_mir.txt", col.names=T, row.names=F, sep="\t", quote=F)
write.table(fltr_rcc, "fltr_rcc.txt", col.names=T, row.names=F, sep="\t", quote=F)
write.table(fltr_sh, "fltr_sh.txt", col.names=T, row.names=F, sep="\t", quote=F)

snp_lq<- unique(select(fltr_lq, snpID, ref, alel))
colnames(snp_lq)[3]="LEIQONG"

snp_mir<- unique(select(fltr_mir, snpID, ref, alel))
colnames(snp_mir)[3]="MIRKADIM"

snp_rcc<- unique(select(fltr_rcc, snpID, ref, alel))
colnames(snp_rcc)[3]="RCC"

snp_sh<- unique(select(fltr_sh, snpID, ref, alel))
colnames(snp_sh)[3]="SAHIWAL"

write.table(snp_lq, "snp_lq.txt", col.names=T, row.names=F, sep="\t", quote=F)
write.table(snp_mir, "snp_mir.txt", col.names=T, row.names=F, sep="\t", quote=F)
write.table(snp_rcc, "snp_rcc.txt", col.names=T, row.names=F, sep="\t", quote=F)
write.table(snp_sh, "snp_sh.txt", col.names=T, row.names=F, sep="\t", quote=F)

##total snp
tsnp.zebu<- Reduce(function(x,y) merge(x, y, by = "snpID", all.x = TRUE, all.y = TRUE),
               list(br, brakmas, ch, hn, kk1, kk2,
               kke, lq, mir, rcc, rs, sh))

snp.info<- Reduce(function(x,y) merge(x, y, by = c("snpID", "ref", "alel"), all.x = TRUE, all.y = TRUE),
              list(br, brakmas, ch, hn, kk1, kk2, kke, lq, mir, rcc, rs, sh))


sum_tsnp<- sqldf("select count(snpID) AS count_alel, snpID from gsnp group by snpID")

#find polymorphic snp
sum_geno<- snp_info %>% 
    group_by(snpID) %>% 
    summarise(alel = paste(alel, collapse="/"))

py<- sum_snp %>%
  filter(str_detect(alel, "/"))

mono <- sum_snp %>%
  filter(!str_detect(alel, "/"))

#replace NC to ref
a <- a %>%
  mutate(BRAHMAN= ifelse(BRAHMAN== "NC", ref,  BRAHMAN))

#
w_tsnp[is.na(w_tsnp)]<- "NC"
w_tsnp$countNC<- rowSums(w_tsnp == "NC")

#split column into new rows
rslt <- zebu2vep %>%
  separate_rows(allele, sep = "/")



  
#snphylo.sh -H HapMap_file -t 8

tsnp<- read.table("RnD/livestock/zebu_tsnp.txt", header=T, sep="\t")

snp.info<- select(tsnp, snpID, new_alel)

snp.info$chr_pos<- snp.info$snpID
setDT(snp.info)
snp.info[, c("a", "b", "pos") := tstrsplit(chr_pos, "_", fixed = TRUE)]
snp.info[, chr := paste(a, b, sep = "_")]

colnames(snp.info)[2]="allele"

snp.info<- select(snp.info, snpID, allele, chr, pos)

a<- merge(snp.info, chr, by="chr")

##hapmap

#double the allele in each row
columns_to_transform <- colnames(a)[-1]  # Exclude the first column (snpID)

# Transform the alleles in the selected columns
a[, columns_to_transform] <- lapply(a[, columns_to_transform], function(x) paste0(x, x))


chr_hapmap<- chr %>%
  filter(new_chr != "X") %>%
  filter(new_chr != "Y")

hapmap29chr<- merge(chr_hapmap, snp.info, by="chr")

hapmap_info$strand<- "NA"
hapmap_info$assembly<- "NA"
hapmap_info$center<- "NA"
hapmap_info$protLSID<- "NA"
hapmap_info$assayLSID<- "NA"
hapmap_info$panel<- "NA"
hapmap_info$QCCode<- "NA"

write.table(hapmap29chr, "hapmap29chr.txt", col.names=T, row.names=F, sep="\t")

write.table(a, "hapmap_alel.txt", col.names=T, row.names=F, sep="\t")

hapmap<- merge(hapmap29chr, a, by="snpID")


#summarize allele
# Function to get unique alleles
sum_alel <- function(row) {
  alel <- unique(unlist(row[2:ncol(row)]))
  paste(sort(alel), collapse = ", ")
}

# Apply the function row-wise
sum_alel_5breeds <- hapmap_5breeds %>%
  rowwise() %>%
  mutate(Allele_Summary = sum_alel(cur_data_all()))

# View summary
sum_alel_5breeds[, c("snpID", "Allele_Summary")]


###
genic_vep<- vep_rslt %>%
  filter(Consequence %in% c("missense_variant", "synonymous_variant", "intron_variant", "5_prime_UTR_variant",
                            "3_prime_UTR_variant"))


