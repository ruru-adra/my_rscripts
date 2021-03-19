fpkm_novel<- fpkm_all %>%
  filter(str_detect(geneID, "^Novel"))

fpkm_tableau2<- filter(fpkm_tableau, Measure_Values >= 2)
fpkm_tableau2<- unique(select(fpkm_tableau2, geneID))

fpkm_novel2<- merge(fpkm_novel, fpkm_tableau2, by="geneID")

fpkm_novel_ann<- nvl_ann %>%
  merge(fpkm_novel2, by="geneID") %>%
  select(geneID, chromosome, KEGG.AnnotInfo, RabGAP.TBC.domain.containing.protein) %>%
  unique()

fpkm_novel_ann2<- filter(fpkm_novel_ann, KEGG.AnnotInfo!= "--")

b<- merge(nvl_ann, rapdb, by="rapdb_desc")

nvl_ann<- merge(fpkm_novel2, nvl_ann, by="geneID")


##unique & difference
nvl_bali<- filter(fpkm_novel, BaliBR != "0")
nvl_ph<- filter(fpkm_novel, PH9BR != "0")
nvl_mrm16<- filter(fpkm_novel, MRM16RR != "0")
nvl_q100<- filter(fpkm_novel, MRQ100RR != "0")
nvl_q76<- filter(fpkm_novel, MRQ76WR != "0")
nvl_mr297<- filter(fpkm_novel, MR297WR != "0")

g.bali<- select(nvl_bali, geneID, BaliBR)
g.ph<- select(nvl_ph, geneID, PH9BR)
g.mrm16<- select(nvl_mrm16, geneID, MRM16RR)
g.q100<- select(nvl_q100, geneID, MRQ100RR)
g.mr297<- select(nvl_mr297, geneID, MR297WR)
g.q76<- select(nvl_q76, geneID, MRQ76WR)

snp_gene<- Reduce(function(x,y) merge(x, y, by = "snp_id", all.x = TRUE, all.y = TRUE),
                  list(a_ns, b_ns))

a<- nvl_gene
a$BaliBR<- as.character(a$BaliBR)
a$PH9BR<- as.character(a$PH9BR)
a$MRM16RR<- as.character(a$MRM16RR)
a$MRQ100RR<- as.character(a$MRQ100RR)
a$MR297WR<- as.character(a$MR297WR)
a$MRQ76WR<- as.character(a$MRQ76WR)

a[is.na(a)] <- "NC"

a$count<- rowSums(a == "NC")

uq_gene<- filter(a, count == "5")
sh_gene<- filter(a, count < '5')

uq_q76<- uq_gene %>%
  filter(MRQ76WR != "NC") %>%
  select(geneID, MRQ76WR)

colnames(kegg_q76)[1]="kegg_term"
colnames(kegg_q76)[3]="keggID"
colnames(kegg_q76)[6]="pvalue"
colnames(kegg_q76)[7]="adj"
colnames(kegg_q76)[8]="spID"

pv_q76<- kegg_q76 %>%
  filter(pvalue < 0.05) %>%
  select(kegg_term, keggID, pvalue, adj, spID)

pv_q76<- pv_q76 %>%
  mutate(spID= strsplit(as.character(spID), ";")) %>%
  unnest(spID)




ls_ftb<- uq_gene %>%
  filter(q76 != "NC") %>%
  select(snpID, q76)

ns_q76<- tsnp %>%
  filter(effect__1=="missense_variant") %>%
  select(snp_name__1, var__1, snp_id__1, chr__1, pos__1, ref__1, alel__1, effect__1, OsID__1)

uqns_q76<- uq_gene %>%
  filter(MRQ76 != "NC") %>%
  select(snpID, MRQ76)


snp_fam_q76<- sqldf("select count (distinct snpID) AS count_snpID, Name from
                      famns_q76 group by Name")
snp_fam_q76<- arrange(snp_fam_q76, desc(count_snpID))

b_fam<- merge(snp_fam_q76, snp_fam_q76_osid, by="Name")

lit_fam_mr297<- merge(famns_mr297, litfun, by="OsID")

vep_q76<- uqns_q76 %>%
  merge(ns_q76, by="snpID") %>%
  select(snpID, chr__1, pos__1, ref__1, alel__1)

a<- sqldf("select count (distinct snpID) AS count_snpID, Name from
                      delfam_q76 group by Name")
a<- arrange(a, desc(count_snpID))

b<- sqldf("select count (distinct OsID) AS count_OsID, Name from
                      delfam_q76 group by Name")
b<- arrange(b, desc(count_OsID))

sum_delfam_q76<- merge(a,b,by="Name")

a<- do.call(rbind, str_split(sp_rslt$subject_id, 'sp'))

sp_go2agbase<- sp_ann %>%
  mutate(goID= strsplit(as.character(goID), ";")) %>%
  unnest(goID) %>%
  filter(goID!="NA") %>%
  select(spID, goID) %>%
  unique()

sp_at<- sp_ann %>%
  filter(Organism=="Arabidopsis thaliana (Mouse-ear cress)") %>%
  mutate(goID= strsplit(as.character(goID), ";")) %>%
  unnest(goID) %>%
  filter(goID!="NA") %>%
  select(spID, goID) %>%
  unique()

sp_os<- sp_ann %>%
  filter(Organism=="Oryza sativa subsp. japonica (Rice)") %>%
  mutate(goID= strsplit(as.character(goID), ";")) %>%
  unnest(goID) %>%
  filter(goID!="NA") %>%
  select(spID, goID) %>%
  unique()

#

lst_pcc_fv_ann<- merge(lst_pcc_fv, ann_funrice_rapdb, by="edge") #all fv, no fltr with fpkm value

pcc_fpkm_valid<- merge(lst_pcc_fv_ann,fpkm_fv_valid, by="node")

write.table(lst_pcc_fv_ann, "~/Desktop/ITEM_OBJ2/coexp/Lampiran_lst_fv_pcc_all_ann.txt", col.names=T, row.names=F, sep="\t", quote=F)
write.table(pcc_fpkm_valid, "~/Desktop/pcc_fpkm_67_fv2cytoscape.txt", col.names=T, row.names=F, sep="\t", quote=F)


degree_s<- select(degree, OsID, Degree, geneName, Symbol, rapdb_desc, pcc, p_value)
degree_s<- arrange(degree_s, desc(Degree))

degree_edge<- degree_s[62:24169, ]
degree_top10<- filter(degree_edge, Degree >= 10)

cluego_keggID<- cluego_kegg %>%
  mutate(ID= strsplit(as.character(ID), ",")) %>%
  unnest(ID)
  #filter(goID!="NA") %>%
  #select(spID, goID) %>%
  #unique()
  
rslt_cluego_keggID<- merge(cluego_keggID, r_cluego, by="ID") #merge id cluego with id uniprot

hub_kegg<- merge(degree_edge, rslt_cluego_keggID, by="OsID") #to get which OsID has kegg + degree
hub_kegg<- arrange(hub_kegg, desc(Degree))

lst_hub_gene<- merge(hub_kegg2, pcc_fpkm_valid, by="edge")

lst_hub_gene2<- filter(lst_hub_gene, GOTerm %in% c("Amino sugar and nucleotide sugar metabolism",
                                                   "Galactose metabolism", "Glycolysis / Gluconeogenesis",
                                                   "Pentose phosphate pathway", "Phenylpropanoid biosynthesis",
                                                   "Phenylpropanoid biosynthesis", "Fructose and mannose metabolism",
                                                   "Phenylalanine metabolism",
                                                   "Plant hormone signal transduction",
                                                   "Starch and sucrose metabolism"))

coexp_hub<- unique(select(lst_hub_gene2, nodeID, node, geneName, edge, Symbol, rapdb_desc, Degree, pcc, p_value))

write.table(coexp_hub, "~/Desktop/coexp_hub.txt", col.names=T, row.names=F, sep="\t", quote=F)

coexp_hub_pv<- filter(coexp_hub, p_value < 0.05)

sp_goID<- sp %>%
  mutate(goID= strsplit(as.character(goID), ";")) %>%
  unnest(goID) %>%
  filter(goID!="NA") %>%
  select(OsID, goID) %>%
  unique()

string_exp<- filter(string, experiments != "NA")
string_coexp<- filter(string, coexpression != "NA")
coexp_exp<- merge(string_exp, string_coexp, by="name")

string2<- do.call(rbind, str_split(string$name, '(pp)'))
string<- cbind(string, string2)

#
string_a
coexp_hubGene

lst_match<- merge(coexp_hubGene, string_a, by=c("node", "edge"))
lst_match2<- merge(coexp_exp, coexp_hubGene, by=c("node", "edge"))
ann_lst_match<- merge(lst_match, string, by=c("node_desc", "edge_desc"))
colnames(string)[15]="node_desc"
colnames(string)[16]="edge_desc"

tf2agbase<- tf2agbase %>%
  mutate(goID= strsplit(as.character(goID), ";")) %>%
  unnest(goID) %>%
  filter(goID!="NA") %>%
  select(OsID, Protein.names, goID) %>%
  unique()

go_cc2<- go_cc %>%
  merge(go_all, by=c("goID", "goTerm")) %>%
  select(GO_Type, goID, goTerm, OsID, total_genes) %>%
  unique()

folate_bp<- merge(folate, go_bp2, by="OsID")
folate_mf<- merge(folate, go_mf2, by="OsID")
folate_cc<- merge(folate, go_cc2, by="OsID")

go_mf2<- filter(go_mf2, goTerm != "molecular_function")

#find folate-tf similar go
mf_tfft<- folate_mf %>%
  merge(go_mf2, by=c("goTerm", "goID")) %>%
  select(goTerm, goID, OsID.y, GO_Type.y, total_genes.y) %>%
  unique()

bp_tfft<- merge(bp_tfft, tf_fam, by="OsID")
mf_tfft<- merge(mf_tfft, tf_fam, by="OsID")
cc_tfft<- merge(cc_tfft, tf_fam, by="OsID")

#merge with nwk grn
colnames(osid_go)[1]="edge_OsID"
grn_go<- merge(osid_go, lst_grn, by="edge_OsID")

#snp-grn folate final
snp_grn<- rsnp %>%
  merge(grn_go, by="edge_OsID") %>%
  select(edge_OsID, edge_geneName, Degree, snpID, ref, alel, chr, pos, snpEffect,
         rapdb_desc, bali, ph, mrm16, q100, mr297, q76) %>%
 unique()


#integrate
coexp_hub_ft<- read.table("hub_gene/hub_ft_bg.txt", header=T, sep="\t", quote="")

colnames(gr_snp)[4]="edge"

#fv hub + tf
osid_hubfv<- unique(select(coexp_hub_fv, edge, rapdb_desc)) #264 genes
snp_gg_coexp_hub_fv<- merge(gr_snp, osid_hubfv, by="edge") #176 snps + 71 genes
osid_tffv<- unique(select(coexp_tf_fv, edge_OsID, edge_geneName)) #155 genes
snp_gg_coexp_tf_fv<- merge(gr_snp, osid_tffv, by="edge_OsID") # 89 snps + 47 genes

#folate hub + tf
osid_hubft<- unique(select(s_hub, edge, edge_geneName)) #204 genes
snp_gg_coexp_hub_ft<- merge(gr_snp, osid_hubft, by="edge")# #210 snps + 90 genes
osid_tfft<- unique(select(coexp_tf_folate, edge_OsID, edge_geneName)) #80 genes
snp_gg_coexp_tf_ft<- merge(gr_snp, osid_tfft, by="edge_OsID") #84 snps + 34 genes

#combine with degree+pcc (from obj 2)
nwk_ft<- merge(coexp_tf_folate, snp_gg_coexp_tf_ft, by=c("edge_OsID","edge_geneName"))
nwk_fv<- merge(coexp_tf_fv, snp_gg_coexp_tf_fv, by=c("edge_OsID","edge_geneName"))

nwk_hubfv<- merge(coexp_hub_fv,snp_gg_coexp_hub_fv, by="edge")
nwk_hubft<- merge(s_hub,snp_gg_coexp_hub_ft, by=c("edge", "edge_geneName"))



nwk_ft<- select(nwk_ft, node_OsID, node_geneName, edge_OsID, edge_geneName, pcc, p_value, Degree, snpID, effect, rapdb_desc, ref, alel, chr, pos, bali, ph, mrm16, q100, mr297, q76)




#folate hub selection
colnames(kobas)[1]="pathway_name"
colnames(kobas)[2]="pathwayID"
colnames(kobas)[7]="uniprotID"

ft_kobas<- kobas %>%
  mutate(uniprotID= strsplit(as.character(uniprotID), ";")) %>%
  unnest(uniprotID) %>%
  select(pathway_name, pathwayID, uniprotID) %>%
  unique()

s_hub<- merge(coexp_hub_ft, ft_kobas, by="edge")
s_hub<- arrange(s_hub, desc(Degree)) #list of folate hub gene

#edit rapdb
rapdb_funrice<- rapdb %>%
  merge(funrice, by="OsID") %>%
  select(OsID, Symbol, rapdb_desc)

notIn<- anti_join(rapdb, rapdb_funrice, by="OsID")

#extract FPKM to create heatmap

fpkm_hubfv<- fpkm %>%
  merge(snp_gg_coexp_hub_fv, by="edge" ) %>%
  select(edge, rapdb_desc.x, BaliBR, PH9BR, MRM16RR, MRQ100RR, MR297WR, MRQ76WR) %>%
  unique()

fpkm_fv<- fpkm %>%
  merge(fv, by="node") %>%
  select(node, geneName, BaliBR, PH9BR, MRM16RR, MRQ100RR, MR297WR, MRQ76WR) %>%
  unique()

fpkm_hubfv$type="hub"
fpkm_fv$type="fv"

fpkm_hubft2tableau<- rbind(fpkm_hubft, fpkm_ft)
fpkm_hubfv2tableau<- rbind(fpkm_hubfv, fpkm_fv)

fpkm_tffv<- fpkm %>%
  merge(coexp_tf_fv, by="edge_OsID" ) %>%
  select(edge_OsID, edge_geneName, BaliBR, PH9BR, MRM16RR, MRQ100RR, MR297WR, MRQ76WR) %>%
  unique()

fpkm_tfft<- fpkm %>%
  merge(coexp_tf_folate, by="edge_OsID" ) %>%
  select(edge_OsID, edge_geneName, BaliBR, PH9BR, MRM16RR, MRQ100RR, MR297WR, MRQ76WR) %>%
  unique()



##score
coexp<- rbind(score_coexp_hubft, score_coexp_tfft)
snp<- rbind(score_snp_hubft, score_snp_tfft)

fpkm<- select(score_fpkm, OsID, type)
coexp<- select(coexp, OsID, type)
snp<- select(snp, OsID, type)


score_ft<- Reduce(function(x,y) merge(x, y, by = "OsID", all.x = TRUE, all.y = TRUE),
                  list(fpkm, coexp, snp))


colnames(score_ft)[2]="fpkm"
colnames(score_ft)[3]="coexp"
colnames(score_ft)[4]="snp"

s_ft<- score_ft
score_ft$fpkm<- as.character(score_ft$fpkm)
score_ft$coexp<- as.character(score_ft$coexp)
score_ft$snp<- as.character(score_ft$snp)


score_ft[is.na(score_ft)] <- "NC"

score_ft$countNC<- rowSums(score_ft == "NC")

p_ft<- filter(s_ft, countNC <= 1)

map_gene<- merge(map_gene, rapdb_funrice, by="OsID")
map_gene2<- merge(map_gene2, rapdb_funrice, by="OsID")
map_folate<- merge(map_folate, rapdb_funrice, by="OsID")


coexp_fv<- unique(select(coexp_hub_fv, node, geneName))
coexp_ft<- unique(select(coexp_hub_ft, node, node_geneName))

snp_fv<- unique(select(snp_gr_fv, OsID, geneName))
snp_ft<- unique(select(snp_gr_ft, OsID, geneName))

snp_ft$type="snp"

#deg analysis
br_gg_fv<- merge(deg_br, gg_fv, by="OsID")
bw_gg_fv<- merge(deg_bw, gg_fv, by="OsID")
rw_gg_fv<- merge(deg_rw, gg_fv, by="OsID")

br_gg_fv$deg="br"
bw_gg_fv$deg="bw"
rw_gg_fv$deg="rw"

br_gg_ft<- merge(deg_br, gg_ft, by="OsID")
bw_gg_ft<- merge(deg_bw, gg_ft, by="OsID")
rw_gg_ft<- merge(deg_rw, gg_ft, by="OsID")

br_gg_ft$deg="br"
bw_gg_ft$deg="bw"
rw_gg_ft$deg="rw"

c<- select(rw_gg_ft, OsID, log2FoldChange, significant, fpkm, coexp, snp, countNC, Symbol, 
                   rapdb_desc, family, type, deg)

deg_gg_fv<- rbind(a,b,c)
deg_gg_ft<- rbind(a,b,c)

#bg folate n fv
br_ft<- merge(deg_br, ft, by="OsID")
bw_ft<- merge(deg_bw, ft, by="OsID")
rw_ft<- merge(deg_rw, ft, by="OsID")

br_ft$deg="br"
bw_ft$deg="bw"
rw_ft$deg="rw"

c<- select(rw_ft, OsID, log2FoldChange, significant, fpkm, coexp, snp, countNC, deg)

write.table(deg_fv, "~/Desktop/ITEM_OBJ3/deg/sig_fv.txt", col.names=T, row.names=F, sep="\t", quote=F)
write.table(deg_ft, "~/Desktop/ITEM_OBJ3/deg/sig_ft", col.names=T, row.names=F, sep="\t", quote=F)


#prepare new table for scoring
deg_rw_up<- deg_gg_ft %>%
  filter(deg=="rw") %>%
  filter(log2FoldChange >= 0)
  deg_rw_up$deg_rw="up"

deg_rw_dwn<- deg_gg_ft %>%
  filter(deg=="rw") %>%
  filter(log2FoldChange < 0)
  deg_rw_dwn$deg_rw="down" 
  
deg_rw<- rbind(deg_rw_up, deg_rw_dwn)

a<- unique(select(deg_br, OsID, deg_br))
b<- unique(select(deg_bw, OsID, deg_bw))
c<- unique(select(deg_rw, OsID, deg_rw))

a_b_c<- Reduce(function(x,y) merge(x, y, by = "OsID", all.x = TRUE, all.y = TRUE),
                  list(a, b, c))

lst_deg_ft<- a_b_c

score_fv<- merge(a, lst_deg_fv, by="OsID")
score_ft<- merge(b, lst_deg_ft, by="OsID")

ann_fv_score<- merge(ann_fv_score, score_fv, by=c("OsID", "fpkm", "coexp", "snp", "type"))
ann_ft_score<- merge(ann_ft_score, score_ft, by=c("OsID", "fpkm", "coexp", "snp", "type"))


msg_type<- node %>%
  merge(msg_type, by="display.name") %>%
  merge(score_ft, by="query.term") %>%
  select(canonical.name,database.identifier,query.term,description,display.name,enhancedLabel.Passthrough
         ,type)
  

a<- merge(string_ft_gg_hubtf_allscore, hub_tf, by="query.term")
msg_queryterm<- anti_join(string_ft_gg_hubtf_allscore, a, by="query.term")

a<- select(a, display.name, coexpression, databases,experiments,textmining,interaction,interspecies,edge,                    
neighborhood,score,selected,shared.interaction,shared.name,canonical.name,database.identifier, 
query.term, description, enhancedLabel.Passthrough, type)

string_ft_gg_hubtf_allscore<- rbind(a, msg_queryterm)

#ppi flavonoid string
id_fv_with_node<- merge(ann_fv, node_fv2, by="query.term")

colnames(node_fv2)[1]="node"
ann_edgefv<- merge(edge_fv, node_fv2, by="node")
a<- merge(ann_edgefv, type_fvhubtf, by="query.term")

nwk<- merge(a, id_fv_with_node, by="query.term")
nwk2<- merge(a, id_fv_with_node, by="display.name")


#string evaluation
fv_ppi<- merge(fv_string, fv_node, by="node")
colnames(fv_ppi)[12]="query.term_node"
colnames(fv_node)[3]="edge"
fv_ppi<- merge(fv_ppi, fv_node, by="edge")
colnames(fv_ppi)[15]="query.term_edge"
colnames(fv_ppi)[10]="database.identifier_node"
colnames(fv_ppi)[13]="database.identifier_edge"
colnames(fv_ppi)[11]="description_node"
colnames(fv_ppi)[14]="description_edge"
fv_ppi<- select(fv_ppi, node, query.term_node, description_node, edge, query.term_edge,
                description_edge, interaction, score, shared.name, coexpression, experiments,
                databases, textmining, interaction, database.identifier_node, 
                database.identifier_edge)

colnames(fv_ppi)[1]="node_string"
colnames(fv_ppi)[2]="node"
colnames(fv_ppi)[4]="edge_string"
colnames(fv_ppi)[5]="edge"

a_hub<- select(fv_gg_hub, node, edge)
a_tf<- select(fv_gg_tf, node_OsID, edge_OsID)
a_tf$type="TF"
a_hub$type="hub"
fv_hubtf<- rbind(a_hub, a_tf)

vldt_fv<- merge(fv_ppi, fv_hubtf, by=c("node", "edge")) #5 matches

  
ft_ppi<- merge(ft_string, ft_node, by="node")
colnames(ft_ppi)[12]="query.term_node"
colnames(ft_node)[3]="edge"
ft_ppi<- merge(ft_ppi, ft_node, by="edge")
colnames(ft_ppi)[22]="query.term_edge"
colnames(ft_ppi)[17]="database.identifier_node"
colnames(ft_ppi)[20]="database.identifier_edge"
colnames(ft_ppi)[18]="description_node"
colnames(ft_ppi)[21]="description_edge"

ft_ppi<- select(ft_ppi, node, query.term_node, description_node, edge, query.term_edge,
                description_edge, interaction, score, shared.name, coexpression, experiments,
                databases, textmining, interaction, database.identifier_node, 
                database.identifier_edge)

colnames(ft_ppi)[1]="node_string"
colnames(ft_ppi)[2]="node"
colnames(ft_ppi)[4]="edge_string"
colnames(ft_ppi)[5]="edge"

vldt_ft<- merge(ft_ppi, ft_hubtf, by=c("node", "edge")) #5 matches



#format table nwk
h_a<- unique(select(ft_gg_hub, node, geneName, edge, Symbol, famili))
h_b<- unique(select(ft_gg_hub, edge, Symbol, node, geneName, famili))
h_a$type="folate"
h_b$type="hub"
h_a$display.name<- h_a$geneName
h_b$display.name<- h_b$Symbol
colnames(h_b)[1]="node"
colnames(h_b)[2]="geneName"
colnames(h_b)[3]="edge"
colnames(h_b)[4]="Symbol"
h_ab<- rbind(h_a, h_b)

tf_a<- unique(select(ft_gg_tf, node, geneName, edge, Symbol, famili))
tf_b<- unique(select(ft_gg_tf, edge, Symbol, node, geneName, famili))
tf_a$type="folate"
tf_b$type="tf"
colnames(tf_b)[1]="node"
colnames(tf_b)[2]="geneName"
colnames(tf_b)[3]="edge"
colnames(tf_b)[4]="Symbol"
tf_a$display.name<- tf_a$geneName
tf_b$display.name<- tf_b$geneName
tf_ab<- rbind(tf_a, tf_b)

nwk_ft_hubtf<- rbind(h_ab, tf_ab)
write.table(nwk_ft_hubtf, "../nwk_ft_hubtf2cys.txt", col.names=T, row.names=F, sep="\t", quote=F)


h_c<- unique(select(nwk_ft_hub_snp, snpID, effect, edge, Symbol, ref, alel, chr, pos, bali, ph, mrm16, q100, mr297, q76))
h_c$type="SNP"
h_c$display.name=h_c$snpID
#h_snp 

tf_c<- unique(select(nwk_ft_tf_snp, snpID, effect, edge, Symbol, ref, alel, chr, pos, bali, ph, mrm16, q100, mr297, q76))
tf_c$type="SNP"
tf_c$display.name<- tf_c$snpID

nwk_ft_hubtf_snp<- unique(rbind(h_c, tf_c))
write.table(nwk_ft_hubtf_snp, "../nwk_ft_gg_hubtf_snp2cys.txt", col.names=T, row.names=F, sep="\t", quote=F)


#snp folate merge
lst_gr_fv<- nwk_fv_all_snp2 %>%
  filter(type=="flavonoid") %>%
  select(node, geneName) %>%
  unique()

lst_gr_fv<- merge(gr_snp, lst_gr_fv, by="OsID")

#merge all hub-tf-snp

h_c<- unique(select(nwk_ft_hub_snp, node, geneName, edge, Symbol, snpID, effect, ref, alel, chr, pos, bali, ph, mrm16, q100, mr297, q76))
h_c$type="folate"
h_c$display.name=h_c$geneName
h_d<- unique(select(nwk_ft_hub_snp, edge, Symbol,node, geneName,snpID, effect, ref, alel, chr, pos, bali, ph, mrm16, q100, mr297, q76))
h_d$type="hub"
h_d$display.name=h_d$Symbol
h_e<- unique(select(nwk_ft_hub_snp, snpID, effect, edge, Symbol, ref, alel, chr, pos, bali, ph, mrm16, q100, mr297, q76))
h_e$type="SNP"
h_e$display.name=h_$snpID









