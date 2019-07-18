## @@ @@ @@ @@ @@ @@ @@ @@ @@ @@  @@ @@ @@ @@ @@ @@ @@ @@ @@ ##
## Generate all the tables for the disease PERCEPTION portal ##
## @@ @@ @@ @@ @@ @@ @@ @@ @@ @@  @@ @@ @@ @@ @@ @@ @@ @@ @@ ##
if ("DB"%in%list.files("./") == FALSE){dir.create("DB")} 
pacientes<-read.csv2("Patients_information_table.txt",stringsAsFactors = F,sep="\t",header=T)
drogas<-readRDS("LINCS_results.rds")
load("Data/Drugs.Rdata") #drugs

# disease_group.tsv
disgroup<-1:length(unique(pacientes$Disease_group_name))
names(disgroup)<-unique(pacientes$Disease_group_name)
disease_group_tsv<-cbind(as.numeric(disgroup),names(disgroup))
colnames(disease_group_tsv)<-c("id","name")
write.table(disease_group_tsv,"DB/disease_group.tsv",quote=F,sep="\t",row.names=F)

# disease_group_properties.tsv
color<-pacientes$Disease_group_color
names(color)<-pacientes$Disease_group_name
disease_group_properties_tsv<-cbind(as.numeric(disgroup),"color",as.character(color[names(disgroup)]))
colnames(disease_group_properties_tsv)<-c("disease_group_id","property","value")
write.table(disease_group_properties_tsv,"DB/disease_group_properties.tsv",quote=F,sep="\t",row.names=F)

# disease.tsv
disease<-1:length(unique(pacientes$Disease_name))
names(disease)<-unique(pacientes$Disease_name)
dis_group<-pacientes$Disease_group_name
names(dis_group)<-pacientes$Disease_name
disease_tsv<-cbind(as.character(disease),names(disease),as.character(disgroup[as.character(dis_group[names(disease)])]))
colnames(disease_tsv)<-c("id","name","disease_group_id")
write.table(disease_tsv,"DB/disease.tsv",quote=F,sep="\t",row.names=F)

# disease_digraph.tsv
disnet<-read.csv2("Networks/Relative_Risks_diseases.txt",stringsAsFactors = F,sep="\t",header=T)
namdis<-pacientes$Disease_name
names(namdis)<-pacientes$Disease
disease_digraph_tsv<-cbind(1:length(disnet[,1]),as.character(disease[as.character(namdis[disnet[,1]])]),
                           as.character(disease[as.character(namdis[disnet[,2]])]),as.numeric(disnet[,3])*as.numeric(disnet[,4]))
colnames(disease_digraph_tsv)<-c("id","disease_a_id","disease_b_id","relative_risk")
write.table(disease_digraph_tsv,"DB/disease_digraph.tsv",quote=F,sep="\t",row.names=F)

# disease_properties.tsv
icd9<-gsub("ICD9_","",pacientes$ICD9)
names(icd9)<-pacientes$Disease_name
icd10<-pacientes$ICD10
names(icd10)<-pacientes$Disease_name
disease_properties_tsv<-cbind(as.character(disease),as.character(icd9[names(disease)]),as.character(icd10[names(disease)]))
colnames(disease_properties_tsv)<-c("disease_id","property","value")
write.table(disease_properties_tsv,"DB/disease_properties.tsv",quote=F,sep="\t",row.names=F)
disdos<-c()
for(a in 1:length(disease)){
  disdos<-rbind(disdos,cbind(as.character(disease[a]),"icd9",as.character(icd9[names(disease[a])])),
                cbind(as.character(disease[a]),"icd10",as.character(icd10[names(disease[a])])),
                cbind(as.character(disease[a]),"color",as.character(color[as.character(dis_group[names(disease[a])])])))
}
colnames(disdos)<-c("disease_id","property","value")
write.table(disdos,"DB/disease_properties.tsv",quote=F,sep="\t",row.names=F)

# patient_subgroup.tsv
patsub<-pacientes$Disease_name
names(patsub)<-pacientes$Cluster_info_cluster_specific
subgrupos<-1:length(unique(pacientes$Cluster_info_cluster_specific))
names(subgrupos)<-unique(pacientes$Cluster_info_cluster_specific)
patient_subgroup_tsv<-cbind(as.character(subgrupos),names(subgrupos),as.numeric(disease[as.character(patsub[names(subgrupos)])]))
colnames(patient_subgroup_tsv)<-c("id","name","disease_id")
write.table(patient_subgroup_tsv,"DB/patient_subgroup.tsv",quote=F,sep="\t",row.names=F)

# patient_subgroup_digraph.tsv
subnet<-read.csv2("Networks/Relative_Risks_subgroups.txt",stringsAsFactors = F,sep="\t",header=T)
patient_subgroup_digraph_tsv<-cbind(1:length(subnet[,1]),as.character(subgrupos[subnet[,1]]),as.character(subgrupos[subnet[,2]]),as.numeric(subnet[,3])*as.numeric(subnet[,4]))
colnames(patient_subgroup_digraph_tsv)<-c("id","patient_subgroup_a_id","patient_subgroup_b_id","relative_risk")
write.table(patient_subgroup_digraph_tsv,"DB/patient_subgroup_digraph.tsv",quote=F,sep="\t",row.names=F)

# study.tsv
splited<-strsplit(pacientes$Patient,"_");studio<-c();for(a in 1:length(splited)){studio<-c(studio,splited[[a]][2])}
study<-1:length(unique(studio))
names(study)<-unique(studio)
study_tsv<-cbind(as.character(study),names(study))
colnames(study_tsv)<-c("id","geo_arrayexpress_code")
write.table(study_tsv,"DB/study.tsv",quote=F,sep="\t",row.names=F)

# patient.tsv
patient_tsv<-cbind(1:length(pacientes$Cluster_info_cluster_specific),as.character(subgrupos[pacientes$Cluster_info_cluster_specific]),as.character(study[studio]))
colnames(patient_tsv)<-c("id","patient_subgroup_id","study_id")
write.table(patient_tsv,"DB/patient.tsv",quote=F,sep="\t",row.names=F)

# patient_graph.tsv
patnet<-read.csv2("Networks/Patient_patient_similarity_network.txt",stringsAsFactors = F,sep="\t",header=T)
pat<-1:length(pacientes$Patient)
names(pat)<-pacientes$Patient
patient_graph_tsv<-cbind(1:length(patnet[,1]),as.character(pat[patnet[,1]]),as.character(pat[patnet[,2]]),patnet[,3])
colnames(patient_graph_tsv)<-c("id","patient_a_id","patient_b_id","interaction_sign")
write.table(patient_graph_tsv,"DB/patient_graph.tsv",quote=F,sep="\t",row.names=F)

# drug.tsv
drogis<-1:length(as.character(unique(drogas[,3])))
names(drogis)<-as.character(unique(drogas[,3]))
drug_tsv<-cbind(as.character(drogis),names(drogis))
colnames(drug_tsv)<-c("id","name")
write.table(drug_tsv,"DB/drug.tsv",quote=F,sep="\t",row.names=F)

# patient_drug_maps.tsv
drojas<-c()
for(a in 1:length(names(drugs))){
  if(length(drugs[[a]]$Similar)>0 && length(drugs[[a]]$Opposite)>0){
    drojas<-rbind(drojas,cbind(names(drugs)[a],drugs[[a]]$Similar,1),cbind(names(drugs)[a],drugs[[a]]$Opposite,-1))
  }
  if(length(drugs[[a]]$Similar)>0 && length(drugs[[a]]$Opposite)==0){
    drojas<-rbind(drojas,cbind(names(drugs)[a],drugs[[a]]$Similar,1))
  }
  if(length(drugs[[a]]$Similar)==0 && length(drugs[[a]]$Opposite)>0){
    drojas<-rbind(drojas,cbind(names(drugs)[a],drugs[[a]]$Opposite,-1))
  }
}
patient_drug_maps_tsv<-cbind(1:length(drojas[,1]),as.character(pat[drojas[,1]]),as.character(drogis[drojas[,2]]),drojas[,3])
colnames(patient_drug_maps_tsv)<-c("id","patient_id","drug_id","regulation_sign")
write.table(patient_drug_maps_tsv,"DB/patient_drug_maps.tsv",quote=F,sep="\t",row.names=F)

# gene.tsv
gene_tsv<-read.csv2("DB/gene.tsv",stringsAsFactors = F,sep="\t",header=T)
# write.table(gene_tsv,"DB/gene.tsv",quote=F,sep="\t",row.names=T)

# patient_gene_maps.tsv
load("Data/Top_500_sDEGs_names.Rdata")
genes<-gene_tsv[,1]
names(genes)<-gene_tsv[,2]

genos<-c()
for(a in 1:length(degs)){genos<-rbind(genos,cbind(names(degs)[a],degs[[a]]$up,1),cbind(names(degs)[a],degs[[a]]$down,-1))}
patient_gene_maps_tsv<-cbind(1:length(genos[,1]),as.character(pat[genos[,1]]),as.character(genes[genos[,2]]),genos[,3])
colnames(patient_gene_maps_tsv)<-c("id","patient_id","gene_id","regulation_sign")
write.table(patient_gene_maps_tsv,"DB/patient_gene_maps.tsv",quote=F,sep="\t",row.names=T)

# patient_subgroup_gene_intersect.tsv
load("Data/Subgroups_genes.Rdata") #commongenes
patient_subgroup_gene_intersect<-c()
for(a in 1:length(names(commongenes))){
  # a<-1
  if(length(commongenes[[a]]$up)>0){
    patient_subgroup_gene_intersect<-rbind(patient_subgroup_gene_intersect,cbind(as.character(subgrupos[names(commongenes)[a]]),as.character(genes[commongenes[[a]]$up]),1))
  }
  if(length(commongenes[[a]]$down)>0){
    patient_subgroup_gene_intersect<-rbind(patient_subgroup_gene_intersect,cbind(as.character(subgrupos[names(commongenes)[a]]),as.character(genes[commongenes[[a]]$down]),-1))
  }
}
patient_subgroup_gene_intersect_tsv<-cbind(1:length(patient_subgroup_gene_intersect[,1]),patient_subgroup_gene_intersect)
colnames(patient_subgroup_gene_intersect_tsv)<-c("id","patient_subgroup_id","gene_id","regulation_sign")
write.table(patient_subgroup_gene_intersect_tsv,"DB/patient_subgroup_gene_intersect.tsv",quote=F,sep="\t",row.names=T)

# patient_subgroup_drug_intersect.tsv
load("Data/Subgroups_drogas.Rdata") #clusterl
patient_subgroup_drug_intersect<-c()
for(a in 1:length(names(clusterl))){
  # a<-1
  if(length(clusterl[[a]]$Up_Dir)>0){
    patient_subgroup_drug_intersect<-rbind(patient_subgroup_drug_intersect,cbind(as.character(subgrupos[names(clusterl)[a]]),as.character(drogis[clusterl[[a]]$Up_Dir]),1))
  }
  if(length(clusterl[[a]]$Down_Inv)>0){
    patient_subgroup_drug_intersect<-rbind(patient_subgroup_drug_intersect,cbind(as.character(subgrupos[names(clusterl)[a]]),as.character(drogis[clusterl[[a]]$Down_Inv]),-1))
  }
}
patient_subgroup_drug_intersect_tsv<-cbind(1:length(patient_subgroup_drug_intersect[,1]),patient_subgroup_drug_intersect)
colnames(patient_subgroup_drug_intersect_tsv)<-c("id","patient_subgroup_id","drug_id","regulation_sign")
write.table(patient_subgroup_drug_intersect_tsv,"DB/patient_subgroup_drug_intersect.tsv",quote=F,sep="\t",row.names=T)

# disease_digraph_properties
## ICD9 ##
icd9<-gsub("ICD9_","",pacientes$ICD9)
names(icd9)<-pacientes$Disease_name
load("Data/PDN_shared_diseases.Rdata") #pdn
uno9<-as.character(icd9[names(disease)[as.numeric(disease_digraph_tsv[,2])]])
dos9<-as.character(icd9[names(disease)[as.numeric(disease_digraph_tsv[,3])]])
mios9<-paste("ICD9-",uno9,"_ICD9-",dos9,sep="")
direc<-as.numeric(disease_digraph_tsv[,4])
conf9<-c() ; for(a in 1:length(mios9)){if(length(which(pdn==mios9[a]))>0 && direc[a]>0){conf9<-c(conf9,a)}}
icd9si<-rep("no",length(disease_digraph_tsv[,1])) ; icd9si[conf9]<-"yes"
## ICD10 ##
icd10<-pacientes$ICD10
names(icd10)<-pacientes$Disease_name
load("Data/DTP_shared_diseases.Rdata") #dtp
uno10<-as.character(icd10[names(disease)[as.numeric(disease_digraph_tsv[,2])]])
dos10<-as.character(icd10[names(disease)[as.numeric(disease_digraph_tsv[,3])]])
mios10<-paste(uno10,"_",dos10,sep="")
direc<-as.numeric(disease_digraph_tsv[,4])
conf10<-c() ; for(a in 1:length(mios10)){if(length(which(dtp==mios10[a]))>0 && direc[a]>0){conf10<-c(conf10,a)}}
icd10si<-rep("no",length(disease_digraph_tsv[,1])) ; icd10si[conf10]<-"yes"
disease_digraph_properties<-cbind(disease_digraph_tsv[,1],icd9si,icd10si)
colnames(disease_digraph_properties)<-c("disease_interaction_id","Hidalgo2009","Jensen2014")
digraph_properties_disease<-c()
appears_in<-c()
for(a in 1:length(disease_digraph_properties[,1])){
  if(disease_digraph_properties[a,2]=="yes" && disease_digraph_properties[a,3]=="yes"){appears_in<-rbind(appears_in,c(a,"appears_in","Hidalgo2009 Jensen2014"))}
  if(disease_digraph_properties[a,2]=="yes" && disease_digraph_properties[a,3]=="no"){appears_in<-rbind(appears_in,c(a,"appears_in","Hidalgo2009"))}
  if(disease_digraph_properties[a,2]=="no" && disease_digraph_properties[a,3]=="yes"){appears_in<-rbind(appears_in,c(a,"appears_in","Jensen2014"))}
}
colnames(appears_in)<-c("disease_digraph_id","property","value")
write.table(appears_in,"DB/disease_digraph_properties.tsv",quote=F,sep="\t",row.names=F)

## supporting_papers
supporting_papers<-rbind(c("Hidalgo2009","19360091","10.1371/journal.pcbi.1000353"),
                         c("Jensen2014","24959948","10.1038/ncomms5022"))
colnames(supporting_papers)<-c("id","pubmed_id","doi")
write.table(supporting_papers,"DB/supporting_papers.tsv",quote=F,sep="\t",row.names=F)


