genic_bali<- filter(c_bali, effect %in% c("missense_variant", "synonymous_variant", "3_prime_UTR_variant",
                                           "5_prime_UTR_variant", "intron_variant"))
                    s_g.bali<- select(genic_bali, chr, pos, ref, alel, snp_id, snp_name)
alel_g.bali<- unique(s_g.bali)

genic_ph<- filter(c_ph, effect %in% c("missense_variant", "synonymous_variant", "3_prime_UTR_variant",
                                      "5_prime_UTR_variant", "intron_variant"))

s_g.ph<- select(genic_ph, chr, pos, ref, alel, snp_id, snp_name)
alel_g.ph<- unique(s_g.ph)

genic_mrm16<- filter(c_mrm16, effect %in% c("missense_variant", "synonymous_variant", "3_prime_UTR_variant",
                                            "5_prime_UTR_variant", "intron_variant"))

s_g.mrm16<- select(genic_mrm16, chr, pos, ref, alel, snp_id, snp_name)
alel_g.mrm16<- unique(s_g.mrm16)

genic_q100<- filter(c_q100, effect %in% c("missense_variant", "synonymous_variant", "3_prime_UTR_variant",
                                          "5_prime_UTR_variant", "intron_variant"))

s_g.q100<- select(genic_q100, chr, pos, ref, alel, snp_id, snp_name)
alel_g.q100<- unique(s_g.q100)

genic_mr297<- filter(c_mr297, effect %in% c("missense_variant", "synonymous_variant", "3_prime_UTR_variant",
                                          "5_prime_UTR_variant", "intron_variant"))
s_g.mr297<- select(genic_mr297, chr, pos, ref, alel, snp_id, snp_name)
alel_g.mr297<- unique(s_g.mr297)

genic_q76<- filter(c_q76, effect %in% c("missense_variant", "synonymous_variant", "3_prime_UTR_variant",
                                            "5_prime_UTR_variant", "intron_variant"))
s_g.q76<- select(genic_q76, chr, pos, ref, alel, snp_id, snp_name)
alel_g.q76<- unique(s_g.q76)

intron_bali<- filter(genic_bali, effect == "intron_variant")
cds_bali<- filter(genic_bali, effect %in% c("missense_variant", "synonymous_variant"))
utr_bali<- filter(genic_bali, effect %in% c("3_prime_UTR_variant", "5_prime_UTR_variant"))
s_i.bali<- select(intron_bali, chr, pos, ref, alel, snp_id, snp_name)
s_cds.bali<- select(cds_bali, chr, pos, ref, alel, snp_id, snp_name)
s_utr.bali<- select(utr_bali, chr, pos, ref, alel, snp_id, snp_name)
alel_i.bali<- unique(s_i.bali)
alel_cds.bali<- unique(s_cds.bali)
alel_utr.bali<- unique(s_utr.bali)

intron_ph<- filter(genic_ph, effect == "intron_variant")
cds_ph<- filter(genic_ph, effect %in% c("missense_variant", "synonymous_variant"))
utr_ph<- filter(genic_ph, effect %in% c("3_prime_UTR_variant", "5_prime_UTR_variant"))
s_i.ph<- select(intron_ph, chr, pos, ref, alel, snp_id, snp_name)
s_cds.ph<- select(cds_ph, chr, pos, ref, alel, snp_id, snp_name)
s_utr.ph<- select(utr_ph, chr, pos, ref, alel, snp_id, snp_name)
alel_i.ph<- unique(s_i.ph)
alel_cds.ph<- unique(s_cds.ph)
alel_utr.ph<- unique(s_utr.ph)

intron_mrm16<- filter(genic_mrm16, effect == "intron_variant")
cds_mrm16<- filter(genic_mrm16, effect %in% c("missense_variant", "synonymous_variant"))
utr_mrm16<- filter(genic_mrm16, effect %in% c("3_prime_UTR_variant", "5_prime_UTR_variant"))
s_i.mrm16<- select(intron_mrm16, chr, pos, ref, alel, snp_id, snp_name)
s_cds.mrm16<- select(cds_mrm16, chr, pos, ref, alel, snp_id, snp_name)
s_utr.mrm16<- select(utr_mrm16, chr, pos, ref, alel, snp_id, snp_name)
alel_i.mrm16<- unique(s_i.mrm16)
alel_cds.mrm16<- unique(s_cds.mrm16)
alel_utr.mrm16<- unique(s_utr.mrm16)

intron_q100<- filter(genic_q100, effect == "intron_variant")
cds_q100<- filter(genic_q100, effect %in% c("missense_variant", "synonymous_variant"))
utr_q100<- filter(genic_q100, effect %in% c("3_prime_UTR_variant", "5_prime_UTR_variant"))
s_i.q100<- select(intron_q100, chr, pos, ref, alel, snp_id, snp_name)
s_cds.q100<- select(cds_q100, chr, pos, ref, alel, snp_id, snp_name)
s_utr.q100<- select(utr_q100, chr, pos, ref, alel, snp_id, snp_name)
alel_i.q100<- unique(s_i.q100)
alel_cds.q100<- unique(s_cds.q100)
alel_utr.q100<- unique(s_utr.q100)

ph_ns<- filter(genic_ph, effect == "missense_variant")
s_ns.ph<- select(ph_ns, chr, pos, ref, alel, snp_id, snp_name)
alel_ns.ph<- unique(s_ns.ph)
ph_syn<- filter(genic_ph, effect == "synonymous_variant")
s_syn.ph<- select(ph_syn, chr, pos, ref, alel, snp_id, snp_name)
alel_syn.ph<- unique(s_syn.ph)

mrm16_ns<- filter(genic_mrm16, effect == "missense_variant")
s_ns.mrm16<- select(mrm16_ns, chr, pos, ref, alel, snp_id, snp_name)
alel_ns.mrm16<- unique(s_ns.mrm16)
mrm16_syn<- filter(genic_mrm16, effect == "synonymous_variant")
s_syn.mrm16<- select(mrm16_syn, chr, pos, ref, alel, snp_id, snp_name)
alel_syn.mrm16<- unique(s_syn.mrm16)

q100_ns<- filter(genic_q100, effect == "missense_variant")
s_ns.q100<- select(q100_ns, chr, pos, ref, alel, snp_id, snp_name)
alel_ns.q100<- unique(s_ns.q100)
q100_syn<- filter(genic_q100, effect == "synonymous_variant")
s_syn.q100<- select(q100_syn, chr, pos, ref, alel, snp_id, snp_name)
alel_syn.q100<- unique(s_syn.q100)

#run total snp

snp_pigment<- Reduce(function(x,y) merge(x, y, by = "snp_id", all.x = TRUE, all.y = TRUE),
                 list(s_snp.bali, s_snp.ph, s_snp.mrm16, s_snp.q100))
snp<- snp_pigment

snp_pigment$bali<- as.character(snp_pigment$bali)
snp_pigment$ph<- as.character(snp_pigment$ph)
snp_pigment$mrm16<- as.character(snp_pigment$mrm16)
snp_pigment$q100<- as.character(snp_pigment$q100)

snp_pigment[is.na(snp_pigment)]<- "nc"
snp_pigment$count<- rowSums(snp_pigment == "nc")
snp_pigment$A<- rowSums(snp_pigment == "A")
snp_pigment$C<- rowSums(snp_pigment == "C")
snp_pigment$G<- rowSums(snp_pigment == "G")
snp_pigment$T<- rowSums(snp_pigment == "T")
snp_pigment$countZero<- rowSums(snp_pigment == "0")
py_pigment<- filter(snp_pigment, countZero == "2")

#create vcf for vep
py_snp<- select(py2vep, snp_id, bali, ph, mrm16, q100)
py_snp$freq<- "2"
freq.py_snp<- py_snp[rep(row.names(py_snp), py_snp$freq),]
tsl_py4<- freq.py_snp[, 1:5]
tsl_py4.ped<- t(tsl_py4)

data<- read.table("data.txt", header=F, sep="\t")
data$pop<- "-9"
data$pop2<- "-9"
data$pop3<- "-9"
data$pop4<- "-9"
data$pop5<- "-9"
colnames(data)[1]="var"
s_data<- select(data, pop, var, pop2, pop3, pop4, pop5)
tsl_v4.bind<- cbind(s_data, tsl_py4.ped)
new<- tsl_v4.bind[2:5 ,]
write.table(new, "tsl_v4.txt", col.names=F, row.names=F, sep="\t", quote=F)

pos_v4<- select(tsl_py4, snp_id)
map_v4<- separate(data=pos_v4 , col=snp_id, into=c("chr", "pos"), sep="\\_")
map_v4<- cbind(pos_v4, map_v4)
map_v4$pop<- '-9'
map_v4$chr<- str_replace_all(map_v4$chr, "chr0", "")
map_v4$chr<- str_replace_all(map_v4$chr, "chr", "")
write.table(map_v4, "tsl_mapV4.txt", col.names=F, row.names=F, sep="\t", quote=F)

#count unique snp
snp_uq<- snp
snp_uq$bali<- as.character(snp_uq$bali)
snp_uq$ph<- as.character(snp_uq$ph)
snp_uq$mrm16<- as.character(snp_uq$mrm16)
snp_uq$q100<- as.character(snp_uq$q100)
snp_uq[is.na(snp_uq)]<- "nc"
snp_uq$count<- rowSums(snp_uq == "nc")

uq_bali<- filter(uq_pigment, bali != "nc")
uq_ph<- filter(uq_pigment, ph != "nc")


#fv
fv_pigment<- Reduce(function(x,y) merge(x, y, by = "snp_id", all.x = TRUE, all.y = TRUE),
                    list(alel_fv.bali, alel_fv.ph, alel_fv.mrm16, alel_fv.q100))

fv_snp<- fv_pigment
fv_snp$bali<- as.character(fv_snp$bali)
fv_snp$ph<- as.character(fv_snp$ph)
fv_snp$mrm16<- as.character(fv_snp$mrm16)
fv_snp$q100<- as.character(fv_snp$q100)
fv_snp[is.na(fv_snp)]<- "nc"
fv_snp$A<- rowSums(fv_snp == "A")
fv_snp$C<- rowSums(fv_snp == "C")
fv_snp$G<- rowSums(fv_snp == "G")
fv_snp$T<- rowSums(fv_snp == "T")
fv_snp$countZero<- rowSums(fv_snp == "0")
py_pigment<- filter(snp_pigment, countZero == "2")
fv_snp$countNC<- rowSums(fv_snp == "nc")






fv_snp[is.na(fv_snp)]<- "nc"
fv_snp$count<- rowSums(fv_snp == "nc")
uq_bali<- filter(uq_pigment, bali != "nc")
uq_ph<- filter(uq_pigment, ph != "nc")

#find ann unique genes
ann_ls.ph<- merge(genic_ph, pv_ph, by ="snp_id")
OsID_alel.q100<- select(ann_ls.q100, OsID, q100)

ann_pigment<- Reduce(function(x,y) merge(x, y, by = "OsID", all.x = TRUE, all.y = TRUE),
                    list(OsID_alel.bali, OsID_alel.ph, OsID_alel.mrm16, OsID_alel.q100))

ann_uq<- ann_pigment
ann_uq$bali<- as.character(ann_uq$bali)
ann_uq$ph<- as.character(ann_uq$ph)
ann_uq$mrm16<- as.character(ann_uq$mrm16)
ann_uq$q100<- as.character(ann_uq$q100)

ann_uq[is.na(ann_uq)]<- "nc"
ann_uq$countNC<- rowSums(ann_uq == "nc")
uq_gene<- filter(ann_uq, countNC == "3") #select unique gene per var
uq_gene.bali<- filter(uq_gene, bali != "nc" )
uq_gene.ph<- filter(uq_gene, ph != "nc" )
uq_gene.mrm16<- filter(uq_gene, mrm16 != "nc" )
uq_gene.q100<- filter(uq_gene, q100 != "nc" )

fv_snp$A<- rowSums(fv_snp == "A")
fv_snp$C<- rowSums(fv_snp == "C")
fv_snp$G<- rowSums(fv_snp == "G")
fv_snp$T<- rowSums(fv_snp == "T")
fv_snp$countZero<- rowSums(fv_snp == "0")
py_pigment<- filter(snp_pigment, countZero == "2")

#go analysis

cluego_q100<- read_excel("../../../../ClueGOResultTable-q100.xls", sheet=2)
bp_q100<- read_excel("../../../../ClueGOResultTable-q100.xls", sheet=4)
mf_q100<- read_excel("../../../../ClueGOResultTable-q100.xls", sheet=5)

colnames(bp_q100)[1]="GOID"
colnames(mf_q100)[1]="GOID"

top_bp.q100<- bp_q100[1:10 ,]
top_mf.q100<- mf_q100[1:10 ,]
merge_topBP.q100<- merge(cluego_q100, top_bp.q100, by="GOID")
merge_topMF.q100<- merge(cluego_q100, top_mf.q100, by="GOID")
top_GO.q100<- rbind(merge_topBP.q100, merge_topMF.q100)
write.table(top_GO.mrm16, "topGO_q100.txt", col.names=T, row.names=F, sep="\t", quote = F)

#friday-13 oct
all_fv<- read.table("fv_genic_all_lv.txt", header=T, sep="\t", quote="")




#annovar
anvr_bali<- read.table("BALI.avinput.variant_function.snp", header=T, sep="\t", fill=TRUE)
colnames(anvr_bali)[4]="pos"
fv_anvr.bali<- merge(anvr_bali, bali_fv, by=c("chr", "pos"))
fv_anvr_genic.bali<- filter(fv_anvr.bali, annoPos %in% c("exonic", "UTR3", "UTR5", "intronic" ))
fv_anvr_genic.bali$anvr_snpName<- paste(fv_anvr_genic.bali$chr, fv_anvr_genic.bali$pos, sep="_")

anvr_ph<- read.table("PH9.avinput.variant_function.snp", header=T, sep="\t", fill=TRUE)
colnames(anvr_ph)[4]="pos"
fv_anvr.ph<- merge(anvr_ph, ph_fv, by=c("chr", "pos"))
fv_anvr_genic.ph<- filter(fv_anvr.ph, annoPos %in% c("exonic", "UTR3", "UTR5", "intronic" ))
fv_anvr_genic.ph$anvr_snpName<- paste(fv_anvr_genic.ph$chr, fv_anvr_genic.ph$pos, sep="_")

anvr_mrm16<- read.table("MRM16.avinput.variant_function.snp", header=T, sep="\t", fill=TRUE)
colnames(anvr_mrm16)[4]="pos"
fv_anvr.mrm16<- merge(anvr_mrm16, mrm16_fv, by=c("chr", "pos"))
fv_anvr_genic.mrm16<- filter(fv_anvr.mrm16, annoPos %in% c("exonic", "UTR3", "UTR5", "intronic" ))
fv_anvr_genic.mrm16$anvr_snpName<- paste(fv_anvr_genic.mrm16$chr, fv_anvr_genic.mrm16$pos, sep="_")

anvr_q100<- read.table("MRQ100.avinput.variant_function.snp", header=T, sep="\t", fill=TRUE)
colnames(anvr_q100)[4]="pos"
fv_anvr.q100<- merge(anvr_q100, q100_fv, by=c("chr", "pos"))
fv_anvr_genic.q100<- filter(fv_anvr.q100, annoPos %in% c("exonic", "UTR3", "UTR5", "intronic" ))
fv_anvr_genic.q100$anvr_snpName<- paste(fv_anvr_genic.q100$chr, fv_anvr_genic.q100$pos, sep="_")

anvr_mr297<- read.table("MR297.avinput.variant_function.snp", header=T, sep="\t", fill=TRUE)
colnames(anvr_mr297)[4]="pos"
fv_anvr.mr297<- merge(anvr_mr297, mr297_fv, by=c("chr", "pos"))
fv_anvr_genic.mr297<- filter(fv_anvr.mr297, annoPos %in% c("exonic", "UTR3", "UTR5", "intronic" ))
fv_anvr_genic.mr297$anvr_snpName<- paste(fv_anvr_genic.mr297$chr, fv_anvr_genic.mr297$pos, sep="_")

anvr_q76<- read.table("Q76.avinput.variant_function.snp", header=T, sep="\t", fill=TRUE)
colnames(anvr_q76)[4]="pos"
fv_anvr.q76<- merge(anvr_q76, q76_fv, by=c("chr", "pos"))
fv_anvr_genic.q76<- filter(fv_anvr.q76, annoPos %in% c("exonic", "UTR3", "UTR5", "intronic" ))
fv_anvr_genic.q76$anvr_snpName<- paste(fv_anvr_genic.q76$chr, fv_anvr_genic.q76$pos, sep="_")

#read new fv
bali_fv<- read.table("bali_genic_fv.txt", header=T, sep="\t", fill=TRUE)
ph_fv<- read.table("ph_genic_fv.txt", header=T, sep="\t", fill=TRUE)
mrm16_fv<- read.table("mrm16_genic_fv.txt", header=T, sep="\t", fill=TRUE)
q100_fv<- read.table("q100_genic_fv.txt", header=T, sep="\t", fill=TRUE)
mr297_fv<- read.table("mr297_genic_fv.txt", header=T, sep="\t", fill=TRUE)
q76_fv<- read.table("q76_genic_fv.txt", header=T, sep="\t", fill=TRUE)

s_alel.bali<- select(ns_fv.bali, snp_id, alel)
s_alel.ph<- select(ns_fv.ph, snp_id, alel)
s_alel.mrm16<- select(ns_fv.mrm16, snp_id, alel)
s_alel.q100<- select(ns_fv.q100, snp_id, alel)

alel_fv_pigment<- Reduce(function(x,y) merge(x, y, by = "snp_id", all.x = TRUE, all.y = TRUE),
                 list(s_alel.bali, s_alel.ph, s_alel.mrm16, s_alel.q100))

alel_fv_pigment$bali<- as.character(alel_fv_pigment$bali)
alel_fv_pigment$ph<- as.character(alel_fv_pigment$ph)
alel_fv_pigment$mrm16<- as.character(alel_fv_pigment$mrm16)
alel_fv_pigment$mr297<- as.character(alel_fv_pigment$mr297)


alel_fv_pigment[is.na(alel_fv_pigment)]<- "nc"
alel_fv_pigment$count<- rowSums(alel_fv_pigment == "nc")
alel_lv$A<- rowSums(alel_lv == "A")
alel_lv$C<- rowSums(alel_lv == "C")
alel_lv$G<- rowSums(alel_lv == "G")
alel_lv$T<- rowSums(alel_lv == "T")