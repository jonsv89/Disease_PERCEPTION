setwd("./")
library("gplots")
## Load the function to extract the genes and drugs involved in comorbidity relations
subgroup_sharing<-function(lista,tabla,cunineg,cunipos,cenf,pathen){
  clusinf<-tabla$Cluster_info_cluster_specific  ;  names(clusinf)<-tabla$Patient  ;  splited<-strsplit(names(clusinf[names(lista)]),"_")
  patname<-c();for(a in 1:length(splited)){patname<-c(patname,paste(as.character(clusinf[names(lista)])[a],splited[[a]][2],splited[[a]][3],splited[[a]][4],sep="_"))}
  names(lista)<-patname
  ## Extract the drugs/genes shared by each subgroup ##
  negativast<-cunineg  ;  positivast<-cunipos
  clusters<-colnames(cenf)  ;  clusterl<-list()
  for(a in 1:length(clusters)){
    coger<-grep(paste("^",clusters[a],"_",sep=""),names(lista))  ;  dirl<-list();invl<-list()
    for(b in 1:length(coger)){drags<-lista[[coger[b]]]
    if(length(grep("down",names(lista[[coger[b]]])))>0){dirl[[b]]<-drags[[2]]  ;  invl[[b]]<-drags[[1]]}
    if(length(grep("Opposite",names(lista[[coger[b]]])))>0){dirl[[b]]<-drags[[1]]  ;  invl[[b]]<-drags[[2]]}
    }
    clusterl[[clusters[a]]]$"Up_Dir"<-Reduce(intersect,dirl)  ;  clusterl[[clusters[a]]]$"Down_Inv"<-Reduce(intersect,invl)
  }
  negativost<-c()
  for(a in 1:length(negativast[,1])){
    unodir<-clusterl[[negativast[a,1]]]$Up_Dir  ;  unoinv<-clusterl[[negativast[a,1]]]$Down_Inv
    dosdir<-clusterl[[negativast[a,2]]]$Up_Dir  ;  dosinv<-clusterl[[negativast[a,2]]]$Down_Inv
    dd<-intersect(unodir,dosdir)  ;  ii<-intersect(unoinv,dosinv)  ;  di<-intersect(unodir,dosinv)  ;  id<-intersect(unoinv,dosdir)
    negativost<-rbind(negativost,c(negativast[a,],length(unodir),length(unoinv),length(dosdir),length(dosinv),length(dd),length(ii),length(di),length(id)))
  }
  colnames(negativost)<-c("Cluster_1","Cluster_2","RR","Direct_1","Inverse_1","Direct_2","Inverse_2","Dir-Dir","Inv-Inv","Dir-Inv","Inv-Dir")
  cumplenunoneg<-intersect(which(as.numeric(negativost[,10])>0),intersect(which(as.numeric(negativost[,8])==0),which(as.numeric(negativost[,9])==0)))
  cumplendosneg<-intersect(which(as.numeric(negativost[,11])>0),intersect(which(as.numeric(negativost[,8])==0),which(as.numeric(negativost[,9])==0)))
  cumplenneg<-unique(c(cumplenunoneg,cumplendosneg))
  intnegt<-negativost[cumplenneg,]
  dn1<-intersect(which(negativost[,5]=="0"),which(negativost[,4]=="0"))
  dn2<-intersect(which(negativost[,7]=="0"),which(negativost[,6]=="0"))
  dsneg<-unique(c(dn1,dn2))
  positivost<-c()
  for(a in 1:length(positivast[,1])){
    unodir<-clusterl[[positivast[a,1]]]$Up_Dir  ;  unoinv<-clusterl[[positivast[a,1]]]$Down_Inv
    dosdir<-clusterl[[positivast[a,2]]]$Up_Dir  ;  dosinv<-clusterl[[positivast[a,2]]]$Down_Inv
    dd<-intersect(unodir,dosdir)  ;  ii<-intersect(unoinv,dosinv)  ;  di<-intersect(unodir,dosinv)  ;  id<-intersect(unoinv,dosdir)
    positivost<-rbind(positivost,c(positivast[a,],length(unodir),length(unoinv),length(dosdir),length(dosinv),length(dd),length(ii),length(di),length(id)))
  }
  colnames(positivost)<-c("Cluster_1","Cluster_2","RR","Direct_1","Inverse_1","Direct_2","Inverse_2","Dir-Dir","Inv-Inv","Dir-Inv","Inv-Dir")
  cumplenunopos<-intersect(which(as.numeric(positivost[,8])>0),intersect(which(as.numeric(positivost[,10])==0),which(as.numeric(positivost[,11])==0)))
  cumplendospos<-intersect(which(as.numeric(positivost[,9])>0),intersect(which(as.numeric(positivost[,10])==0),which(as.numeric(positivost[,11])==0)))
  cumplenpos<-unique(c(cumplenunopos,cumplendospos))
  intpost<-positivost[cumplenpos,]
  dp1<-intersect(which(positivost[,5]=="0"),which(positivost[,4]=="0"))
  dp2<-intersect(which(positivost[,7]=="0"),which(positivost[,6]=="0"))
  dspos<-unique(c(dp1,dp2))
  ## Represent again the heatmap but in this case only with the clusters sharing drugs in the correct direction
  interestmapt<-matrix(ncol=length(cenf[,1]),nrow=length(cenf[,1]),0)
  rownames(interestmapt)<-rownames(cenf)  ;  colnames(interestmapt)<-colnames(cenf)
  for(a in 1:length(intnegt[,1])){interestmapt[as.character(intnegt[a,1]),as.character(intnegt[a,2])]<-as.numeric(intnegt[a,3])*(-1)}
  for(a in 1:length(intpost[,1])){interestmapt[as.character(intpost[a,1]),as.character(intpost[a,2])]<-as.numeric(intpost[a,3])}
  ## Plotting heatmaps ##
  palette.breaks=c(seq(min(interestmapt),-1,length=100),seq(-0.9999999,0.9999999,length=100),seq(1,max(interestmapt),length=100))
  color.palette<-colorRampPalette(c("#FF0000","#FFFFFF","#1022ea"))(length(palette.breaks)-1)
  splited<-strsplit(colnames(interestmapt),".",fixed=T);indef<-c();num_clus<-c();for(a in 1:length(splited)){indef<-c(indef,splited[[a]][1]);num_clus<-c(num_clus,splited[[a]][2])}
  nombre<-paste(as.character(figuras$Nombres[indef]),"Cluster",num_clus,sep=" ")
  colnames(cenf)<-nombre  ;  rownames(cenf)<-nombre  ;  etiquetas<-as.character(figuras$Colores[as.character(figuras$Enfermedades[indef])])
  interestmapt2<-interestmapt  ;  colnames(interestmapt2)<-colnames(cenf2)  ;  rownames(interestmapt2)<-rownames(cenf2)
  
  rencilla<-rbind(cbind(intnegt[,1:2],-1),cbind(intpost[,1:2],1))
  colnames(rencilla)<-c("Interactor_1","Interactor_2","Interaction")
  red<-rbind(cbind(intpost[,c(1,2,3)],1),cbind(intnegt[,c(1,2,3)],-1))
  colnames(red)<-c("Interactor_1","Interactor_2","Intensity","Direction")
  write.table(red,pathen,row.names=F,sep="\t",quote=T)
  resultados<-list("Clusterl"=clusterl,"Interestmapt2"=interestmapt2,"ColorPalette"=color.palette,"PaletteBreaks"=palette.breaks,"Etiquetas"=etiquetas,"Interestmapt"=interestmapt,
                   "Intnegt"=intnegt,"Intpost"=intpost)
  return(resultados)
}
if ("HeatMaps"%in%list.files("Plots/") == FALSE){dir.create("Plots/HeatMaps")}
if ("Supplementary_documents"%in%list.files("./") == FALSE){dir.create("./Supplementary_documents")}

## @ @ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@@@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @ @ ##
## Association of the number of genes by subgroup: real vs expected by disease ##
## @ @ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@@@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @ @ ##
## This part of the code will tell us which are the subgroups presenting more shared genes than the expected by chance (calculated randomly assigning patients from the same disease
## to a subgroups with the same size of the one we are comparing to)
load("Data/Subgroups_common_genes_shuffled_intra_disease.Rdata") # listorras
load("Data/Subgroups_genes.Rdata") # commongenes
length(names(commongenes))  ## number of subgroups obtained 
info<-read.csv2("Patients_information_table.txt",stringsAsFactors = F,sep="\t",header=T)
culur<-info$Disease_group_color  ;  names(culur)<-info$Disease
diseases<-names(listorras[[1]])
## represent the shuffled ones in boxplots
plotbox<-c()
reals<-c()
for(a in 1:length(diseases)){
  ads<-c()
  for(z in 1:length(names(listorras))){
    ads<-c(ads,length(listorras[[z]][[which(names(listorras[[z]])==diseases[a])]]$up)+length(listorras[[z]][[which(names(listorras[[z]])==diseases[a])]]$down))
  }
  plotbox<-cbind(plotbox,ads)
  ## adding real values
  reals<-c(reals,length(commongenes[[which(names(commongenes)==diseases[a])]]$up)+length(commongenes[[which(names(commongenes)==diseases[a])]]$down))
}
colnames(plotbox)<-diseases
plotbox<-as.data.frame(plotbox)
elmax<-c();for(a in 1:length(plotbox[1,])){elmax<-c(elmax,max(plotbox[,a]))}
## calculate the pvalues for each disease
pvals<-c()
for(a in 1:length(plotbox[1,])){pvals<-c(pvals,length(which(plotbox[,a]>=reals[a]))/length(plotbox[,a]))}
color<-rep("blue",length(reals))
color[which(pvals<=0.05)]<-"red"
maxims<-c();for(a in 1:length(plotbox[1,])){maxims<-c(maxims,max(plotbox[,a]))}
color[setdiff(which(maxims==0),which(pvals<=0.05))]<-"green"
## plot the obtained results
colaxis<-as.character(culur[gsub("\\..+","",colnames(plotbox))])
boxplot(plotbox,pch=20,border=colaxis,las=2,cex=0.1,par(mar=c(15,5,4,2)+0.6,cex.axis=0.1))
text(reals,"*",col=color,cex=0.5)
## Calculate the mean number of subgroups per disease and the mean number of patients per subgroup
averagesize<-table(info[-c(grep("NS",info$Cluster_number_cluster_specific),grep("Cluster_of_size",info$Cluster_info_cluster_specific)),8])
mean(as.numeric(averagesize))  ## mean number of patients per subgroup
mean(as.numeric(table(gsub("\\..+","",names(averagesize)))))  ## mean number of subgroups per disease

                       ## @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ ##
                       ## Genes shared by all the patients within a subgroup in Alzheimer's disease and lung cancer -->  enrichment analyses on them ##
                       ## @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ ##

load("Data/Top_500_sDEGs_names.Rdata") # degs
tabla<-read.csv2("Patients_information_table.txt",stringsAsFactors = F,sep="\t",header=T)
groups<-unique(tabla$Cluster_info_cluster_specific)
if(length(grep("NS$",groups))>0){groups<-groups[-c(grep("NS$",groups),grep("Cluster_of_size",groups))]}  ;  commongenes<-list()
for(a in 1:length(groups)){
  ichoose<-tabla$Patient[which(tabla$Cluster_info_cluster_specific==groups[a])]
  ups<-list();downs<-list()  ;  for(b in 1:length(ichoose)){ups[[ichoose[b]]]<-degs[[ichoose[b]]]$up  ;  downs[[ichoose[b]]]<-degs[[ichoose[b]]]$down}
  commongenes[[groups[a]]]$up<-Reduce(intersect,ups)  ;  commongenes[[groups[a]]]$down<-Reduce(intersect,downs)
}
diseases<-c("Alzheimer","NSCLC")
info<-read.csv2("Patients_information_table.txt",stringsAsFactors = F,sep="\t",header=T)  ;  size<-table(info$Cluster_info_cluster_specific)
subl<-list()  ;  tam<-c()
for(a in 1:length(diseases)){took<-grep(diseases[a],names(commongenes))
for(b in took){subl[[names(commongenes)[b]]]<-commongenes[[b]]  ;  tam<-rbind(tam,c(length(commongenes[[b]]$up),length(commongenes[[b]]$down)))}
}
## Select more genes (there are too few in common): Load biogrid network, select the first neighbour of each gene deregulated in all the patients within each subgroup ##
## @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@  @@ @@ @@ @@ @@ @@ @ @ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ ##
biogrid<-read.csv2("Data/BIOGRID-ALL-3.4.164.tab2.txt",stringsAsFactors = F,sep="\t")
if(length(biogrid[1,])==24){biogrid<-biogrid[intersect(which(biogrid$Organism.Interactor.A==9606),which(biogrid$Organism.Interactor.B==9606)),c(8,9,10,11,12,13)]}
load("Data/ExpressionGenes_Together/Acne_GSE53795.Rdata")  ;  genes<-rownames(expression)
genuko<-strsplit(unique(c(biogrid$Synonyms.Interactor.A,biogrid$Synonyms.Interactor.B)),"|",fixed=T)
kk<-c()  ;  for(a in 1:length(genuko)){ku<-genuko[[a]]  ;  for(b in 1:length(ku)){kk<-c(kk,ku[b])}}
length(intersect(toupper(unique(kk)),toupper(genes)))
comunes<-intersect(toupper(genes),toupper(unique(c(biogrid$Official.Symbol.Interactor.A,biogrid$Official.Symbol.Interactor.B))))
bigrid<-cbind(toupper(biogrid$Official.Symbol.Interactor.A),toupper(biogrid$Official.Symbol.Interactor.B))
#Generate a list of first interactors based on biogrid network
genl<-list()
for(a in 1:length(comunes)){
  ba<-c();if(length(which(bigrid[,1]==comunes[a]))>0){ba<-bigrid[which(bigrid[,1]==comunes[a]),2]}
  be<-c();if(length(which(bigrid[,2]==comunes[a]))>0){be<-bigrid[which(bigrid[,2]==comunes[a]),1]}
  genl[[comunes[a]]]<-unique(c(ba,be))
}
save(genl,file="Data/Biogrid_gene_interactors.Rdata")

## Load the first neighbours of each gene based on biogrid network
load("Data/Biogrid_gene_interactors.Rdata")

## Generate a list with all the subgroups' genes and their first neighbours for enrichment analyses
subl2<-list()
for(a in 1:length(names(subl))){
  upi<-c()  ;  downi<-c()
  if(length(subl[[a]]$up)>0){mis<-subl[[a]]$up  ;  ups<-c()  ;  for(b in 1:length(mis)){ups<-c(ups,genl[[mis[b]]])}  ;  upi<-unique(c(ups,mis))  ;  subl2[[names(subl)[a]]]$up<-upi}
  if(length(upi)==0){subl2[[names(subl)[a]]]$up<-""}
  if(length(subl[[a]]$down)>0){mis<-subl[[a]]$down  ;  downs<-c();for(b in 1:length(mis)){downs<-c(downs,genl[[mis[b]]])}
  downi<-unique(c(downs,mis))  ;  subl2[[names(subl)[a]]]$down<-downi}  ;  if(length(downi)==0){subl2[[names(subl)[a]]]$down<-""}
}
save(subl2,file="Data/Genes_shared_by_subgroups.Rdata")

#Run enrichment analyses
library("gProfileR")
load("Data/Genes_shared_by_subgroups.Rdata")
background<-read.csv2("Background.txt",stringsAsFactors = F,sep="\t",header=F)[,1]
pathways<-list()
for(a in 1:length(subl2)){
  if(length(subl2[[a]]$down)>1){
    pat<-gprofiler(subl2[[a]]$down,organism="hsapiens",custom_bg=background)  ;  pathways[[names(subl2)[a]]]$Down<-pat
    print(paste(names(subl2)[a],length(pat$term.name),"down",sep="  "))
  }
  if(length(subl2[[a]]$up)>1){
    pat<-gprofiler(subl2[[a]]$up,organism="hsapiens",custom_bg=background)  ;  pathways[[names(subl2)[a]]]$Up<-pat
    print(paste(names(subl2)[a],length(pat$term.name),"up",sep="  "))
  }
}
save(pathways,file="Data/Subgroups_pathways.Rdata")

## Load pathways associated to each patient-subgroup
load("Data/Subgroups_pathways.Rdata") #pathways

enri_reac<-list()  ;  enri_iden<-list()  ;  tabla<-c()  ;  topreac<-list()
for(a in 1:length(pathways)){
  dep<-names(pathways[[a]])
  for(b in 1:length(dep)){
    patss<-pathways[[a]][[dep[b]]][,c(9,10,12,14,3)]  ;  interpatss<-patss[c(which(patss[,2]=="rea")),]
    if(sum(dim(interpatss))>5){
      #We introduce a filtering step where we are going to select only those pathways with at least one of our subgroups' genes
      genes<-subl[[a]][[tolower(dep[b])]]  ;  vects<-c()
      for(c in 1:length(genes)){vect<-c();for(d in 1:length(interpatss[,1])){vect<-c(vect,length(which(strsplit(interpatss[d,4],",")[[1]]==genes[c])))};vects<-cbind(vects,vect)}
      suma<-apply(vects,1,sum)
      reacs<-suma[order(suma,decreasing=T)]/length(genes)  ;  names(reacs)<-interpatss[order(suma,decreasing=T),3]
      dreacs<-suma[order(suma,decreasing=T)]/length(genes)  ;  names(dreacs)<-interpatss[order(suma,decreasing=T),1]
      if(dep[b]=="Down"){
        print(paste(names(pathways)[a],"Down",length(genes),length(patss[,1]),length(reacs),length(which(reacs>0)),sep="  "))
        tabla<-rbind(tabla,c(names(pathways)[a],"Down",length(genes),length(patss[,1]),length(reacs),length(which(reacs>0))))}
      if(dep[b]=="Up"){
        print(paste(names(pathways)[a],"Up",length(genes),length(patss[,1]),length(reacs),length(which(reacs>0)),sep="  "))
        tabla<-rbind(tabla,c(names(pathways)[a],"Up",length(genes),length(patss[,1]),length(reacs),length(which(reacs>0))))}
      reacs<-interpatss[order(as.numeric(interpatss[,5]),decreasing=F),5]  ;  names(reacs)<-interpatss[order(as.numeric(interpatss[,5]),decreasing=F),3]
      dreacs<-interpatss[order(as.numeric(interpatss[,5]),decreasing=F),5]  ;  names(dreacs)<-interpatss[order(as.numeric(interpatss[,5]),decreasing=F),1]
      if(dep[b]=="Down"){enri_reac[[names(pathways)[a]]]$Down<-reacs  ;  enri_iden[[names(pathways)[a]]]$Down<-dreacs
      if(length(reacs)<10){topreac[[names(pathways)[a]]]$Down<-reacs}
      if(length(reacs)>=10){topreac[[names(pathways)[a]]]$Down<-reacs[1:10]}
      }
      if(dep[b]=="Up"){enri_reac[[names(pathways)[a]]]$Up<-reacs  ;  enri_iden[[names(pathways)[a]]]$Up<-dreacs}
    }
    if(sum(dim(interpatss))==5){print(paste(names(pathways)[a],dep[b],"No reactome pathways",sep="  "))}
  }
}

## @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @ @ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ ##
## Generate a list with the drugs generating gene-expression changes similar and opposite to the ones observed in our patients ##
## @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @ @ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ ##
tabla<-readRDS("LINCS_results.rds")  ;  patdrugs<-tabla[,-c(1:8)]  ;  drugs<-list()
for(a in 1:length(patdrugs[1,])){
  pos<-unique(as.character(tabla[which(as.numeric(patdrugs[,a])>=2),3]))  ;  neg<-unique(as.character(tabla[which(as.numeric(patdrugs[,a])<=-2),3]))
  similar<-setdiff(pos,neg)  ;  opposite<-setdiff(neg,pos)  ;  drugs[[colnames(patdrugs)[a]]]$Similar<-similar  ;  drugs[[colnames(patdrugs)[a]]]$Opposite<-opposite
}
save(drugs,file="Data/Drugs.Rdata")

## Information for coloring properly the labels
tabla<-read.csv2("Patients_information_table.txt",stringsAsFactors = F,sep="\t",header=T)
nombres<-unique(tabla$Disease_name)  ;  names(nombres)<-unique(tabla$Disease)
color<-unique(tabla$Disease_group_color)
names(color)<-c("Color_12","Color_13","Color_2","Color_7","Color_0","Color_9","Color_6","Color_10","Color_8","Color_5","Color_4","Color_1","Color_14","Color_3","Color_16")
ic<-tabla[,4:3]  ;  ic<-ic[-which(duplicated(ic)),]  ;  icd9s<-ic[,1]  ;  names(icd9s)<-ic[,2]
coli<-names(color)  ;  names(coli)<-as.character(color)  ;  culi<-cbind(tabla[,3],as.character(coli[tabla$Disease_group_color]))  ;  culi<-culi[-which(duplicated(culi)),]
enfermedades<-culi[,2]  ;  names(enfermedades)<-culi[,1]
figuras<-list("Nombres"=nombres,"Colores"=color,"ICD9s"=icd9s,"Enfermedades"=enfermedades)
izena<-tabla$Disease_name  ;  names(izena)<-tabla$Disease

length(grep("Cluster",tabla$Cluster_info_cluster_specific)) ## Patients clustered alone 

## Number of clusters composed by patients from different studies 
tablita<-tabla[-c(grep("NS$",tabla$Cluster_info_cluster_specific),grep("^^",tabla$Cluster_info_cluster_specific,fixed=T)),]
subgroups<-unique(tablita$Cluster_info_cluster_specific)
howmany<-c()
for(a in 1:length(subgroups)){
  splited<-strsplit(tablita$Patient[which(tablita$Cluster_info_cluster_specific==subgroups[a])],"_")
  study<-c();for(b in 1:length(splited)){study<-c(study,splited[[b]][2])}  ;  howmany<-c(howmany,length(unique(study)))
}
length(unique(gsub("\\..+","",subgroups[which(howmany!=1)]))) ## number of diseases composed by subgroups composed by patients from different studies

splited<-strsplit(tabla$Patient,"_");stus<-c();for(a in 1:length(splited)){stus<-c(stus,splited[[a]][2])};stuts<-unique(paste(tabla$Disease,stus,sep="_"))
length(unique(gsub("\\..+","",subgroups[which(howmany!=1)])))/length(which(as.numeric(as.character(table(gsub("_.+","",stuts))))>1))

length(which(howmany>1))/length(howmany)  ## Percentage of true subgroups composed by patients from different studies  --> 20.74176%


## @@ @@ @@ ##
## Diseases ##
## @@ @@ @@ ##
load("Data/RR/Disease_pRR.Rdata")
load("Data/RR/Disease_nRR.Rdata")
RRspositivas<-RRpos[,c(1:2,5)]
RRsnegativas<-RRneg[,c(1:2,5)]
enfermedade<-unique(tabla$Disease[order(tabla$Disease_group_number,decreasing = F)])
## Select the significant ones
rrpossig<-RRspositivas[which(as.numeric(RRspositivas[,3])>=1),]
rrnegsig<-RRsnegativas[which(as.numeric(RRsnegativas[,3])>=1),]
## Remove the overlapping interactions (both up and down)
pospaste<-paste(rrpossig[,1],rrpossig[,2],sep="_") ## 3.090
negpaste<-paste(rrnegsig[,1],rrnegsig[,2],sep="_") ## 2.504
comunes<-intersect(pospaste,negpaste) ## 422
## Generate the unique interactions
indpos<-c()  ;  indneg<-c()  ;  for(a in 1:length(comunes)){indpos<-c(indpos,which(pospaste==comunes[a]))  ;  indneg<-c(indneg,which(negpaste==comunes[a]))}
cunineg<-rrnegsig[-indneg,]  ;  cunipos<-rrpossig[-indpos,]
network<-rbind(cbind(cunipos,1),cbind(cunineg,-1))  ;  colnames(network)<-c("Disease_1","Disease_2","RR","Interaction")
write.table(network,"Networks/Relative_Risks_diseases.txt",quote=F,sep="\t",row.names=F)
length(cunipos[,1])/length(network[,1]) # percentage of positive interactions
group<-tabla$Disease_group_name  ;  names(group)<-tabla$Disease
length(which(as.character(group[cunipos[,1]])==as.character(group[cunipos[,2]])))/length(cunipos[,1])
## Generate the matrix
para_heatpam<-enfermedade  ;  cenf<-matrix(ncol=length(para_heatpam),nrow=length(para_heatpam),0)  ;  colnames(cenf)<-para_heatpam  ;  rownames(cenf)<-para_heatpam
for(a in 1:length(cunipos[,1])){cenf[as.character(cunipos[a,1]),as.character(cunipos[a,2])]<-as.numeric(cunipos[a,3])}
for(a in 1:length(cunineg[,1])){cenf[as.character(cunineg[a,1]),as.character(cunineg[a,2])]<-(-1*as.numeric(cunineg[a,3]))}
## Plotting heatmaps ##
palette.breaks=c(seq(min(cenf),-1,length=100),seq(-0.9999999,0.9999999,length=100),seq(1,max(cenf),length=100))
color.palette<-colorRampPalette(c("#FF0000","#FFFFFF","#1022ea"))(length(palette.breaks)-1)
cenf2<-cenf
enf<-cenf2  ;  save(enf,file="Data/Disease_matrix.Rdata")
etiquetas<-as.character(figuras$Colores[as.character(figuras$Enfermedades[colnames(cenf2)])])
colnames(cenf2)<-as.character(izena[colnames(cenf2)])  ;  rownames(cenf2)<-as.character(izena[rownames(cenf2)])
dev.off()
heatmap.2(cenf2,scale="none",dendrogram = "none",trace="none",col=color.palette,breaks=palette.breaks,
          Colv=FALSE,Rowv=FALSE,cexRow = 0.5,cexCol = 0.5,"key"=FALSE,margins = c(10,15),
          lmat=rbind(c(0,3,0),c(2,1,0),c(0,4,0)),lhei=c(0.5,8,0),lwid=c(0.5,8,0),colRow =etiquetas,colCol = etiquetas)
# Save it as a pdf
pdf(file="Plots/HeatMaps/Diseases.pdf")
heatmap.2(cenf2,scale="none",dendrogram = "none",trace="none",col=color.palette,breaks=palette.breaks,
          Colv=FALSE,Rowv=FALSE,cexRow = 0.2,cexCol = 0.2,"key"=FALSE,margins = c(10,15),
          lmat=rbind(c(0,3,0),c(2,1,0),c(0,4,0)),lhei=c(0.5,8,0),lwid=c(0.5,8,0),colRow =etiquetas,colCol = etiquetas)
dev.off()


## @@ @@ @@ @@ @@ @@ ##
## Patient subgroups ##
## @@ @@ @@ @@ @@ @@ ##
size<-4
## Extract the number of patients per cluster
mayores<-names(table(tabla$Cluster_info_cluster_specific))[which(as.numeric(table(tabla$Cluster_info_cluster_specific))>=size)]  ;  mayores<-mayores[-grep("NS$",mayores)]
names(mayores)<-mayores
load("Data/RR/Patient_subgroup_pRR.Rdata")
load("Data/RR/Patient_subgroup_nRR.Rdata")
RRspositivas<-RRpos[-unique(c(grep("Cluster-of",RRpos[,1]),grep("Cluster-of",RRpos[,2]))),c(1:2,5)]
RRspositivas<-RRspositivas[-unique(c(grep("NS$",RRspositivas[,1]),grep("NS$",RRspositivas[,2]))),]
RRsnegativas<-RRneg[-unique(c(grep("Cluster-of",RRneg[,1]),grep("Cluster-of",RRneg[,2]))),c(1:2,5)]
RRsnegativas<-RRsnegativas[-unique(c(grep("NS$",RRsnegativas[,1]),grep("NS$",RRsnegativas[,2]))),]
enf<-unique(tabla$Cluster_info_cluster_specific[order(tabla$Disease_group_number,decreasing = F)])  ;  enfermedade<-enf[-grep("Cluster",enf)]
## Select the significant ones
rrpossig<-RRspositivas[which(as.numeric(RRspositivas[,3])>=1),]
rrnegsig<-RRsnegativas[which(as.numeric(RRsnegativas[,3])>=1),]
## Remove the overlapping interactions (both up and down)
pospaste<-paste(rrpossig[,1],rrpossig[,2],sep="_") ## 84.059
negpaste<-paste(rrnegsig[,1],rrnegsig[,2],sep="_") ## 24.820
comunes<-intersect(pospaste,negpaste) ## 1.023
## Generate the unique interactions
indpos<-c()  ;  indneg<-c()  ;  for(a in 1:length(comunes)){indpos<-c(indpos,which(pospaste==comunes[a]))  ;  indneg<-c(indneg,which(negpaste==comunes[a]))}
cunineg<-rrnegsig[-indneg,]  ;  cunipos<-rrpossig[-indpos,]
network<-rbind(cbind(cunipos,1),cbind(cunineg,-1))  ;  colnames(network)<-c("Subgroup_1","Subgroup_2","RR","Interaction")
if(size==2){write.table(network,"Networks/Relative_Risks_subgroups.txt",quote=F,sep="\t",row.names=F)}
## Generate the matrix
para_heatpam<-as.character(mayores[enfermedade])
if(length(which(is.na(para_heatpam)))>0){para_heatpam<-para_heatpam[-which(is.na(para_heatpam))]}
cenf<-matrix(ncol=length(para_heatpam),nrow=length(para_heatpam),0)  ;  colnames(cenf)<-para_heatpam  ;  rownames(cenf)<-para_heatpam
## Select only the clusters with the desired size
uno<-which(is.na(as.character(mayores[cunipos[,1]])))  ;  dos<-which(is.na(as.character(mayores[cunipos[,2]])))
if((length(uno)>0 && length(dos)>0)==TRUE){cunipos<-cunipos[-unique(c(uno,dos)),]}
uno<-which(is.na(as.character(mayores[cunineg[,1]])))  ;  dos<-which(is.na(as.character(mayores[cunineg[,2]])))
if((length(uno)>0 && length(dos)>0)==TRUE){cunineg<-cunineg[-unique(c(uno,dos)),]}
for(a in 1:length(cunipos[,1])){cenf[as.character(cunipos[a,1]),as.character(cunipos[a,2])]<-as.numeric(cunipos[a,3])}
for(a in 1:length(cunineg[,1])){cenf[as.character(cunineg[a,1]),as.character(cunineg[a,2])]<-(-1*as.numeric(cunineg[a,3]))}
## Plotting
palette.breaks=c(seq(min(cenf),-1,length=100),seq(-0.9999999,0.9999999,length=100),seq(1,max(cenf),length=100))
color.palette<-colorRampPalette(c("#FF0000","#FFFFFF","#1022ea"))(length(palette.breaks)-1)
cenf2<-cenf
splited<-strsplit(colnames(cenf2),".",fixed=T);indef<-c();num_clus<-c();for(a in 1:length(splited)){indef<-c(indef,splited[[a]][1]);num_clus<-c(num_clus,splited[[a]][2])}
nombre<-paste(as.character(figuras$Nombres[indef]),"Cluster",num_clus,sep=" ")
colnames(cenf2)<-nombre  ;  rownames(cenf2)<-nombre  ;  etiquetas<-as.character(figuras$Colores[as.character(figuras$Enfermedades[indef])])
dev.off()
heatmap.2(cenf2,scale="none",dendrogram = "none",trace="none",col=color.palette,breaks=palette.breaks,
          Colv=FALSE,Rowv=FALSE,cexRow = 0.1,cexCol = 0.1,"key"=FALSE,margins = c(10,15),
          lmat=rbind(c(0,3,0),c(2,1,0),c(0,4,0)),lhei=c(0.5,8,0),lwid=c(0.5,8,0),colRow =etiquetas,colCol = etiquetas)
# Save it as a pdf
pdf(file=paste("Plots/HeatMaps/Cluster_of_size_",size,"_cluster_focused.pdf",sep=""))
heatmap.2(cenf2,scale="none",dendrogram = "none",trace="none",col=color.palette,breaks=palette.breaks,
          Colv=FALSE,Rowv=FALSE,cexRow = 0.1,cexCol = 0.1,"key"=FALSE,margins = c(10,15),
          lmat=rbind(c(0,3,0),c(2,1,0),c(0,4,0)),lhei=c(0.5,8,0),lwid=c(0.5,8,0),colRow =etiquetas,colCol = etiquetas)
dev.off()

## Drugs ##
## @@ @@ ##
load("Data/Drugs.Rdata")
tabla<-read.csv2("Patients_information_table.txt",stringsAsFactors = F,sep="\t",header=T)
## Run the function
clusterl_drugs<-subgroup_sharing(drugs,tabla,cunineg,cunipos,cenf,paste("Networks/Cluster_of_size_",size,"_with_shared_drugs_cluster_focused.txt",sep=""))

if(size==2){clusterl<-clusterl_drugs$Clusterl;save(clusterl,file="Data/Subgroups_drogas.Rdata")}
dev.off()
heatmap.2(clusterl_drugs$Interestmapt2,scale="none",dendrogram = "none",trace="none",col=clusterl_drugs$ColorPalette,breaks=clusterl_drugs$PaletteBreaks,
          Colv=FALSE,Rowv=FALSE,cexRow = 0.1,cexCol = 0.1,"key"=FALSE,margins = c(10,15),
          lmat=rbind(c(0,3,0),c(2,1,0),c(0,4,0)),lhei=c(0.5,8,0),lwid=c(0.5,8,0),colRow =clusterl_drugs$Etiquetas,colCol = clusterl_drugs$Etiquetas)
# Save it as a pdf
pdf(file=paste("Plots/HeatMaps/Cluster_of_size_",size,"_with_shared_drugs_cluster_focused.pdf",sep=""))
heatmap.2(clusterl_drugs$Interestmapt2,scale="none",dendrogram = "none",trace="none",col=clusterl_drugs$ColorPalette,breaks=clusterl_drugs$PaletteBreaks,
          Colv=FALSE,Rowv=FALSE,cexRow = 0.1,cexCol = 0.1,"key"=FALSE,margins = c(10,15),
          lmat=rbind(c(0,3,0),c(2,1,0),c(0,4,0)),lhei=c(0.5,8,0),lwid=c(0.5,8,0),colRow =clusterl_drugs$Etiquetas,colCol = clusterl_drugs$Etiquetas)
dev.off()

## Genes ##
## @@ @@ ##
load("Data/Top_500_sDEGs_names.Rdata")
tabla<-read.csv2("Patients_information_table.txt",stringsAsFactors = F,sep="\t",header=T)
## Run the function
clusterl_genes<-subgroup_sharing(degs,tabla,cunineg,cunipos,cenf,paste("Networks/Cluster_of_size_",size,"_with_shared_genes_cluster_focused.txt",sep=""))
## when run for size 2 clusters ##
# save(clusterl_genes,file="./Supplementary_documents/Clusters_with_associated_genes.Rdata")
save(clusterl_genes,file="./Data/Clusters_by_shared_genes.Rdata")

dev.off()
heatmap.2(clusterl_genes$Interestmapt2,scale="none",dendrogram = "none",trace="none",col=clusterl_genes$ColorPalette,breaks=clusterl_genes$PaletteBreaks,
          Colv=FALSE,Rowv=FALSE,cexRow = 0.1,cexCol = 0.1,"key"=FALSE,margins = c(10,15),
          lmat=rbind(c(0,3,0),c(2,1,0),c(0,4,0)),lhei=c(0.5,8,0),lwid=c(0.5,8,0),colRow =clusterl_genes$Etiquetas,colCol = clusterl_genes$Etiquetas)
# Save it as a pdf
pdf(file=paste("Plots/HeatMaps/Cluster_of_size_",size,"_with_shared_genes_cluster_focused.pdf",sep=""))
heatmap.2(clusterl_genes$Interestmapt2,scale="none",dendrogram = "none",trace="none",col=clusterl_genes$ColorPalette,breaks=clusterl_genes$PaletteBreaks,
          Colv=FALSE,Rowv=FALSE,cexRow = 0.1,cexCol = 0.1,"key"=FALSE,margins = c(10,15),
          lmat=rbind(c(0,3,0),c(2,1,0),c(0,4,0)),lhei=c(0.5,8,0),lwid=c(0.5,8,0),colRow =clusterl_genes$Etiquetas,colCol = clusterl_genes$Etiquetas)
dev.off()
inpos<-clusterl_genes$Intpost[which(gsub("\\..+","",clusterl_genes$Intpost[,1])!=gsub("\\..+","",clusterl_genes$Intpost[,2])),]
inneg<-clusterl_genes$Intnegt[which(gsub("\\..+","",clusterl_genes$Intnegt[,1])!=gsub("\\..+","",clusterl_genes$Intnegt[,2])),]
nut<-rbind(cbind(inpos[,1:3],1),cbind(inneg[,1:3],-1))

length(unique(c(nut[,1],nut[,2])))  ## number of subgroups
length(nut[,1])  ## number of interactions with shared genes
disinter<-unique(paste(gsub("\\..+","",nut[,1]),gsub("\\..+","",nut[,2]),nut[,4],sep="_"))
length(disinter) ## number of interactions between diseases in the size 4 scn
splited<-strsplit(disinter,"_");mandis<-c();for(a in 1:length(splited)){mandis<-c(mandis,splited[[a]][1],splited[[a]][2])}
length(unique(mandis))  ## number of diseases represented by the size 4 scn

disnet<-read.csv2("Networks/Relative_Risks_diseases.txt",stringsAsFactors = F,sep="\t",header=T)
pdisnet<-paste(disnet[,1],disnet[,2],disnet[,4],sep="_")
length(setdiff(disinter,pdisnet))/length(disinter)  ## percentage of newly detected interactions in the size 4 subgroups scn at the disease level not detected in the DMSN
newlydetect<-setdiff(disinter,pdisnet)

## Look for genes involved in newly described subgroup relations ##
## @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ ##
## Positive relations
casos<-clusterl_genes$Intpost[intersect(grep("Asthma",clusterl_genes$Intpost[,1]),grep("NSCLC",clusterl_genes$Intpost[,2])),]
if(sum(dim(casos))>0){
  for(a in 1:length(casos[,1])){
    u<-intersect(clusterl_genes$Clusterl[[which(names(clusterl_genes$Clusterl)==casos[a,1])]]$Up_Dir,clusterl_genes$Clusterl[[which(names(clusterl_genes$Clusterl)==casos[a,2])]]$Up_Dir)
    d<-intersect(clusterl_genes$Clusterl[[which(names(clusterl_genes$Clusterl)==casos[a,1])]]$Down_Inv,clusterl_genes$Clusterl[[which(names(clusterl_genes$Clusterl)==casos[a,2])]]$Down_Inv)
    print(paste(casos[a,1],casos[a,2],sep="-"))
    print(u)
    print(d)
  }
}
if(sum(dim(casos))==0){
  u<-intersect(clusterl_genes$Clusterl[[which(names(clusterl_genes$Clusterl)==casos[1])]]$Up_Dir,clusterl_genes$Clusterl[[which(names(clusterl_genes$Clusterl)==casos[2])]]$Up_Dir)
  d<-intersect(clusterl_genes$Clusterl[[which(names(clusterl_genes$Clusterl)==casos[1])]]$Down_Inv,clusterl_genes$Clusterl[[which(names(clusterl_genes$Clusterl)==casos[2])]]$Down_Inv)
  print(paste(casos[1],casos[2],sep="-"))
  print(u)
  print(d)
}

                                                        ## @@ @@ @@ @@ @@  @@ @@ @@ @@ @@ ##
                                                        ## Check the overlap with the PDN ##
                                                        ## @@ @@ @@ @@ @@  @@ @@ @@ @@ @@ ##
icd9s<-tabla$ICD9 ; names(icd9s)<-tabla$Disease
# pos
plus<-newlydetect[grep("_1",newlydetect)]
pRR<-c();splited<-strsplit(plus,"_");for(a in 1:length(splited)){pRR<-rbind(pRR,c(splited[[a]][1],splited[[a]][2],splited[[a]][3]))}
prr<-cbind(gsub("_","-",as.character(icd9s[pRR[,1]])),gsub("_","-",as.character(icd9s[pRR[,2]])),pRR[,3])
# neg
res<-newlydetect[grep("_-1",newlydetect)]
nRR<-c();splited<-strsplit(res,"_");for(a in 1:length(splited)){nRR<-rbind(nRR,c(splited[[a]][1],splited[[a]][2],splited[[a]][3]))}
nrr<-cbind(gsub("_","-",as.character(icd9s[nRR[,1]])),gsub("_","-",as.character(icd9s[nRR[,2]])),nRR[,3])
husig<-read.csv2("Data/PDN_hidalgo.net",stringsAsFactors = F,sep="\t",header=F)   
network2<-cbind(paste("ICD9-",husig[,1],sep=""),paste("ICD9-",husig[,2],sep=""),husig[,7])
network<-rbind(cbind(as.character(network2[,1]),as.character(network2[,2]),network2[,3]),
               cbind(as.character(network2[,2]),as.character(network2[,1]),network2[,3]))
# Select only those codes analyzed in both diseases
uno<-as.character(network[,1]);dos<-as.character(network[,2])
interactores<-intersect(unique(c(prr[,1],prr[,2],nrr[,1],nrr[,2])),unique(c(uno,dos)))
indu<-c()   ;   indd<-c();for(a in 1:length(interactores)){indu<-c(indu,which(uno==interactores[a]))   ;   indd<-c(indd,which(dos==interactores[a]))}
network2<-cbind(uno[intersect(indu,indd)],dos[intersect(indu,indd)],network[intersect(indu,indd),3])
network2<-network2[which(as.numeric(network2[,3])>1),]
notconnected<-setdiff(interactores,unique(c(network2[,1],network2[,2])))
# Select significant positive and negative RR interactions, removing those in common
posit<-prr[,1:2]   ;   negat<-nrr[,1:2]
indu<-c()   ;   indd<-c();for(a in 1:length(interactores)){indu<-c(indu,which(posit[,1]==interactores[a]))   ;   indd<-c(indd,which(posit[,2]==interactores[a]))}
positivos<-posit[intersect(indu,indd),]
indu<-c()   ;   indd<-c();for(a in 1:length(interactores)){indu<-c(indu,which(negat[,1]==interactores[a]))   ;   indd<-c(indd,which(negat[,2]==interactores[a]))}
negativos<-negat[intersect(indu,indd),]
posipeg<-paste(positivos[,1],positivos[,2],sep="_")   ;   negapeg<-paste(negativos[,1],negativos[,2],sep="_")
netwpeg<-paste(network2[,1],network2[,2],sep="_")
length(posipeg)  # newly detected pRR interactions between diseases
length(intersect(netwpeg,posipeg))  # number of those interactions described by the PDN
length(intersect(netwpeg,posipeg))/length(posipeg) # percentage of newly described interactions previously described by the PDN


                                                       ## @@ @@ @@ @@ @@  @@ @@ @@ @@ @@ ##
                                                       ## Check the overlap with the PDN ##
                                                       ## @@ @@ @@ @@ @@  @@ @@ @@ @@ @@ ##

icd10s<-tabla$ICD10 ; names(icd10s)<-tabla$Disease
# pos
plus<-newlydetect[grep("_1",newlydetect)]
pRR<-c();splited<-strsplit(plus,"_");for(a in 1:length(splited)){pRR<-rbind(pRR,c(splited[[a]][1],splited[[a]][2],splited[[a]][3]))}
prr<-cbind(gsub("_","-",as.character(icd10s[pRR[,1]])),gsub("_","-",as.character(icd10s[pRR[,2]])),pRR[,3])
# neg
res<-newlydetect[grep("_-1",newlydetect)]
nRR<-c();splited<-strsplit(res,"_");for(a in 1:length(splited)){nRR<-rbind(nRR,c(splited[[a]][1],splited[[a]][2],splited[[a]][3]))}
nrr<-cbind(gsub("_","-",as.character(icd10s[nRR[,1]])),gsub("_","-",as.character(icd10s[nRR[,2]])),nRR[,3])

network<-read.csv2("Epidemiological_networks/Disease_pairs.csv",stringsAsFactors = F,sep=",",header=T)[,1:3]

# Select only those codes analyzed in both diseases
uno<-as.character(network[,1]);dos<-as.character(network[,2])
interactores<-intersect(unique(c(prr[,1],prr[,2],nrr[,1],nrr[,2])),unique(c(uno,dos)))
indu<-c()   ;   indd<-c();for(a in 1:length(interactores)){indu<-c(indu,which(uno==interactores[a]))   ;   indd<-c(indd,which(dos==interactores[a]))}
network2<-cbind(uno[intersect(indu,indd)],dos[intersect(indu,indd)],network[intersect(indu,indd),3])
notconnected<-setdiff(interactores,unique(c(network2[,1],network2[,2])))
# Select significant positive and negative RR interactions, removing those in common
posit<-prr[,1:2]   ;   negat<-nrr[,1:2]
indu<-c()   ;   indd<-c();for(a in 1:length(interactores)){indu<-c(indu,which(posit[,1]==interactores[a]))   ;   indd<-c(indd,which(posit[,2]==interactores[a]))}
positivos<-posit[intersect(indu,indd),]
indu<-c()   ;   indd<-c();for(a in 1:length(interactores)){indu<-c(indu,which(negat[,1]==interactores[a]))   ;   indd<-c(indd,which(negat[,2]==interactores[a]))}
negativos<-negat[intersect(indu,indd),]
posipeg<-paste(positivos[,1],positivos[,2],sep="_")   ;   negapeg<-paste(negativos[,1],negativos[,2],sep="_")
netwpeg<-paste(network2[,1],network2[,2],sep="_")
length(posipeg)  # newly detected pRR interactions between diseases
length(intersect(netwpeg,posipeg))  # number of those interactions described by the PDN
length(intersect(netwpeg,posipeg))/length(posipeg) # percentage of newly described interactions previously described by the PDN

                                                                  ## @@ @@ @@ @@ @@ @@ @@ @ @@ @@ @@ @@ @@ @@ @@ @@ ##
                                                                  ## Generate heatmaps for the diseases of interest ##
                                                                  ## @@ @@ @@ @@ @@ @@ @@ @ @@ @@ @@ @@ @@ @@ @@ @@ ##

## Color the labels properly
colorines<-tabla$Disease_group_color
names(colorines)<-paste(tabla$Disease_name," Cluster ",gsub(".+\\.","",tabla$Cluster_info_cluster_specific),sep="")
## Disease of interest
diseases<-c("Alzheimer","NSCLC")
## Just diseases ##
## @@ @@ @ @@ @@ ##
load("Data/Disease_matrix.Rdata")
indcol<-c()  ;  indrow<-c()  ;  for(a in 1:length(diseases)){indcol<-c(indcol,grep(diseases[a],colnames(enf)))  ;  indrow<-c(indrow,grep(diseases[a],rownames(enf)))}
displo<-enf[indrow,indcol]  ;  colnames(displo)<-as.character(izena[colnames(displo)])  ;  rownames(displo)<-as.character(izena[rownames(displo)])
dev.off()
heatmap.2(displo,scale="none",dendrogram = "none",trace="none",col=color.palette,breaks=palette.breaks,
          Colv=FALSE,Rowv=FALSE,cexRow = 0.8,cexCol = 0.8,"key"=FALSE,margins = c(10,15),
          lmat=rbind(c(0,3,0),c(2,1,0),c(0,4,0)),lhei=c(0.5,8,0),lwid=c(0.5,8,0))
pdf("Plots/HeatMaps/AD_NSCLC_disease.pdf")
heatmap.2(displo,scale="none",dendrogram = "none",trace="none",col=color.palette,breaks=palette.breaks,
          Colv=FALSE,Rowv=FALSE,cexRow = 0.8,cexCol = 0.8,"key"=FALSE,margins = c(10,15),
          lmat=rbind(c(0,3,0),c(2,1,0),c(0,4,0)),lhei=c(0.5,8,0),lwid=c(0.5,8,0))
dev.off()

## Just clusters ##
## @ @ @@ @@ @ @ ##
# diseases<-c("Alzheimer","Schiz")
indcol<-c()  ;  indrow<-c()  ;  for(a in 1:length(diseases)){indcol<-c(indcol,grep(diseases[a],colnames(cenf2)))  ;  indrow<-c(indrow,grep(diseases[a],rownames(cenf2)))}
etiquetas<-as.character(colorines[colnames(cenf2[indrow,indcol])])
dev.off()
heatmap.2(cenf2[indrow,indcol],scale="none",dendrogram = "none",trace="none",col=color.palette,breaks=palette.breaks,
          Colv=FALSE,Rowv=FALSE,cexRow = 0.5,cexCol = 0.5,"key"=FALSE,margins = c(10,15),
          lmat=rbind(c(0,3,0),c(2,1,0),c(0,4,0)),lhei=c(0.5,8,0),lwid=c(0.5,8,0),colRow =etiquetas,colCol = etiquetas)
pdf("Plots/HeatMaps/AD_NSCLC_subgroups_cluster_focused.pdf")
heatmap.2(cenf2[indrow,indcol],scale="none",dendrogram = "none",trace="none",col=color.palette,breaks=palette.breaks,
          Colv=FALSE,Rowv=FALSE,cexRow = 0.5,cexCol = 0.5,"key"=FALSE,margins = c(10,15),
          lmat=rbind(c(0,3,0),c(2,1,0),c(0,4,0)),lhei=c(0.5,8,0),lwid=c(0.5,8,0),colRow =etiquetas,colCol = etiquetas)
dev.off()

## Shared Drugs ##
## @@ @@  @@ @@ ##
indcol<-c()  ;  indrow<-c()
for(a in 1:length(diseases)){indcol<-c(indcol,grep(diseases[a],colnames(clusterl_drugs$Interestmapt2))) ; indrow<-c(indrow,grep(diseases[a],rownames(clusterl_drugs$Interestmapt2)))}
etiquetas<-as.character(colorines[colnames(clusterl_drugs$Interestmapt2[indrow,indcol])])
dev.off()
heatmap.2(clusterl_drugs$Interestmapt2[indrow,indcol],scale="none",dendrogram = "none",trace="none",col=clusterl_drugs$ColorPalette,breaks=clusterl_drugs$PaletteBreaks,
          Colv=FALSE,Rowv=FALSE,cexRow = 0.5,cexCol = 0.5,"key"=FALSE,margins = c(10,15),
          lmat=rbind(c(0,3,0),c(2,1,0),c(0,4,0)),lhei=c(0.5,8,0),lwid=c(0.5,8,0),colRow =etiquetas,colCol = etiquetas)
pdf("Plots/HeatMaps/AD_NSCLC_with_shared_drugs_cluster_focused.pdf")
heatmap.2(clusterl_drugs$Interestmapt2[indrow,indcol],scale="none",dendrogram = "none",trace="none",col=clusterl_drugs$ColorPalette,breaks=clusterl_drugs$PaletteBreaks,
          Colv=FALSE,Rowv=FALSE,cexRow = 0.5,cexCol = 0.5,"key"=FALSE,margins = c(10,15),
          lmat=rbind(c(0,3,0),c(2,1,0),c(0,4,0)),lhei=c(0.5,8,0),lwid=c(0.5,8,0),colRow =etiquetas,colCol = etiquetas)
dev.off()

## Shared Genes ##
## @@ @@  @@ @@ ##
indcol<-c()  ;  indrow<-c()
for(a in 1:length(diseases)){indcol<-c(indcol,grep(diseases[a],colnames(clusterl_genes$Interestmapt2))) ; indrow<-c(indrow,grep(diseases[a],rownames(clusterl_genes$Interestmapt2)))}
colorines<-tabla$Disease_group_color
names(colorines)<-paste(tabla$Disease_name," Cluster ",gsub(".+\\.","",tabla$Cluster_info_cluster_specific),sep="")
etiquetas<-as.character(colorines[colnames(clusterl_genes$Interestmapt2[indrow,indcol])])
dev.off()
heatmap.2(clusterl_genes$Interestmapt2[indrow,indcol],scale="none",dendrogram = "none",trace="none",col=clusterl_genes$ColorPalette,breaks=clusterl_genes$PaletteBreaks,
          Colv=FALSE,Rowv=FALSE,cexRow = 0.5,cexCol = 0.5,"key"=FALSE,margins = c(10,15),
          lmat=rbind(c(0,3,0),c(2,1,0),c(0,4,0)),lhei=c(0.5,8,0),lwid=c(0.5,8,0),colRow =etiquetas,colCol = etiquetas)
pdf("Plots/HeatMaps/AD_NSCLC_with_shared_genes_cluster_focused.pdf")
heatmap.2(clusterl_genes$Interestmapt2[indrow,indcol],scale="none",dendrogram = "none",trace="none",col=clusterl_genes$ColorPalette,breaks=clusterl_genes$PaletteBreaks,
          Colv=FALSE,Rowv=FALSE,cexRow = 0.5,cexCol = 0.5,"key"=FALSE,margins = c(10,15),
          lmat=rbind(c(0,3,0),c(2,1,0),c(0,4,0)),lhei=c(0.5,8,0),lwid=c(0.5,8,0),colRow =etiquetas,colCol = etiquetas)
dev.off()

                                  ## @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ ##
                                  ##Extract the genes involved in the interactions of interest for pathways analyses ##
                                  ## @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ ##


con1<-c()  ;  con2<-c()  ;  cop1<-c()  ;  cop2<-c()
for(a in 1:length(diseases)){
  con1<-c(con1,grep(diseases[a],clusterl_genes$Intnegt[,1]))  ;  con2<-c(con2,grep(diseases[a],clusterl_genes$Intnegt[,2]))
  cop1<-c(cop1,grep(diseases[a],clusterl_genes$Intpost[,1]))  ;  cop2<-c(cop2,grep(diseases[a],clusterl_genes$Intpost[,2]))
}
save(clusterl_genes,file="Data/Clusters_by_shared_genes.Rdata")
intneg<-clusterl_genes$Intnegt[intersect(con1,con2),]  ;  intpos<-clusterl_genes$Intpost[intersect(cop1,cop2),]
# intneg<-clusterl_genes$Intnegt  ;  intpos<-clusterl_genes$Intpost
## Extract the genes involved in the interactions of interest
if(length(which(gsub("\\..+","",intneg[,1])==gsub("\\..+","",intneg[,2])))>0){intneg2<-intneg[-which(gsub("\\..+","",intneg[,1])==gsub("\\..+","",intneg[,2])),]}
if(length(which(gsub("\\..+","",intneg[,1])==gsub("\\..+","",intneg[,2])))==0){intneg2<-intneg}
sharpats<-list()
for(a in 1:length(intneg2[,1])){
  genis<-c(intersect(clusterl_genes$Clusterl[[which(names(clusterl_genes$Clusterl)==intneg2[a,1])]]$Up_Dir,
                     clusterl_genes$Clusterl[[which(names(clusterl_genes$Clusterl)==intneg2[a,2])]]$Down_Inv),
           intersect(clusterl_genes$Clusterl[[which(names(clusterl_genes$Clusterl)==intneg2[a,1])]]$Down_Inv,
                     clusterl_genes$Clusterl[[which(names(clusterl_genes$Clusterl)==intneg2[a,2])]]$Up_Dir))
  print(paste(intneg2[a,1]," ",intneg2[a,2]," -> ",genis,sep=""))
  
  updown<-intersect(clusterl_genes$Clusterl[[which(names(clusterl_genes$Clusterl)==intneg2[a,1])]]$Up_Dir,
                    clusterl_genes$Clusterl[[which(names(clusterl_genes$Clusterl)==intneg2[a,2])]]$Down_Inv)
  downup<-intersect(clusterl_genes$Clusterl[[which(names(clusterl_genes$Clusterl)==intneg2[a,1])]]$Down_Inv,
                    clusterl_genes$Clusterl[[which(names(clusterl_genes$Clusterl)==intneg2[a,2])]]$Up_Dir)
  sharpats[[paste(intneg2[a,1],intneg2[a,2],sep="_")]]$Updown<-updown
  sharpats[[paste(intneg2[a,1],intneg2[a,2],sep="_")]]$Downup<-downup
}
if(length(which(gsub("\\..+","",intpos[,1])==gsub("\\..+","",intpos[,2])))>0){intpos2<-intpos[-which(gsub("\\..+","",intpos[,1])==gsub("\\..+","",intpos[,2])),]}
if(length(which(gsub("\\..+","",intpos[,1])==gsub("\\..+","",intpos[,2])))==0){intpos2<-intpos}

if(length(intpos2==11)){
  genis<-c(intersect(clusterl_genes$Clusterl[[which(names(clusterl_genes$Clusterl)==intpos2[1])]]$Up_Dir,
                     clusterl_genes$Clusterl[[which(names(clusterl_genes$Clusterl)==intpos2[2])]]$Up_Dir),
           intersect(clusterl_genes$Clusterl[[which(names(clusterl_genes$Clusterl)==intpos2[1])]]$Down_Inv,
                     clusterl_genes$Clusterl[[which(names(clusterl_genes$Clusterl)==intpos2[2])]]$Down_Inv))
  print(paste(intpos2[1]," ",intpos2[2]," -> ",genis,sep=""))
  upup<-intersect(clusterl_genes$Clusterl[[which(names(clusterl_genes$Clusterl)==intpos2[1])]]$Up_Dir,
                  clusterl_genes$Clusterl[[which(names(clusterl_genes$Clusterl)==intpos2[2])]]$Up_Dir)
  downdown<-intersect(clusterl_genes$Clusterl[[which(names(clusterl_genes$Clusterl)==intpos2[1])]]$Down_Inv,
                      clusterl_genes$Clusterl[[which(names(clusterl_genes$Clusterl)==intpos2[2])]]$Down_Inv)
  sharpats[[paste(intpos2[1],intpos2[2],sep="_")]]$UpUp<-upup
  sharpats[[paste(intpos2[1],intpos2[2],sep="_")]]$DownDown<-downdown
}
for(a in 1:length(intpos2[,1])){
  genis<-c(intersect(clusterl_genes$Clusterl[[which(names(clusterl_genes$Clusterl)==intpos2[a,1])]]$Up_Dir,
                     clusterl_genes$Clusterl[[which(names(clusterl_genes$Clusterl)==intpos2[a,2])]]$Up_Dir),
           intersect(clusterl_genes$Clusterl[[which(names(clusterl_genes$Clusterl)==intpos2[a,1])]]$Down_Inv,
                     clusterl_genes$Clusterl[[which(names(clusterl_genes$Clusterl)==intpos2[a,2])]]$Down_Inv))
  print(paste(intpos2[a,1]," ",intpos2[a,2]," -> ",genis,sep=""))
  upup<-intersect(clusterl_genes$Clusterl[[which(names(clusterl_genes$Clusterl)==intpos2[a,1])]]$Up_Dir,
                    clusterl_genes$Clusterl[[which(names(clusterl_genes$Clusterl)==intpos2[a,2])]]$Up_Dir)
  downdown<-intersect(clusterl_genes$Clusterl[[which(names(clusterl_genes$Clusterl)==intpos2[a,1])]]$Down_Inv,
                    clusterl_genes$Clusterl[[which(names(clusterl_genes$Clusterl)==intpos2[a,2])]]$Down_Inv)
  sharpats[[paste(intpos2[a,1],intpos2[a,2],sep="_")]]$UpUp<-upup
  sharpats[[paste(intpos2[a,1],intpos2[a,2],sep="_")]]$DownDown<-downdown
}
save(sharpats,file="Data/AD_NSCLC_overlapping_genes.Rdata")
load("Data/AD_NSCLC_overlapping_genes.Rdata") #sharpats

## Select more genes (there are too few in common): Load biogrid network, select the first neighbour of each gene deregulated in all the patients within each subgroup ##
## @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@  @@ @@ @@ @@ @@ @@ @ @ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ ##
load("Data/Biogrid_gene_interactors.Rdata")

#Generate a list with all the subgroups' genes and their first neighbours for enrichment analyses
sharpats2<-list()
for(a in 1:length(names(sharpats))){
  upi<-c()  ;  downi<-c()
  if(length(sharpats[[a]]$Updown)>0){mis<-sharpats[[a]]$Updown  ;  ups<-c()  ;  for(b in 1:length(mis)){ups<-c(ups,genl[[mis[b]]])}
  upi<-unique(c(ups,mis))  ;  sharpats2[[names(sharpats)[a]]]$Updown<-upi}
  if(length(upi)==0){sharpats2[[names(sharpats)[a]]]$Updown<-""}
  if(length(sharpats[[a]]$Downup)>0){mis<-sharpats[[a]]$Downup  ;  downs<-c();for(b in 1:length(mis)){downs<-c(downs,genl[[mis[b]]])}
  downi<-unique(c(downs,mis))  ;  sharpats2[[names(sharpats)[a]]]$Downup<-downi}  ;  if(length(downi)==0){sharpats2[[names(sharpats)[a]]]$Downup<-""}
  
  ## Relacion directa
  if(length(sharpats[[a]]$DownDown)>0){mis<-sharpats[[a]]$DownDown  ;  ups<-c()  ;  for(b in 1:length(mis)){ups<-c(ups,genl[[mis[b]]])}
  upi<-unique(c(ups,mis))  ;  sharpats2[[names(sharpats)[a]]]$Updown<-upi}
}

save(sharpats2,file="Data/AD_NSCLC_overlapping_genes.Rdata")
load("Data/AD_NSCLC_overlapping_genes.Rdata") #sharpats2
library("gProfileR")
background<-read.csv2("Data/Background.txt",stringsAsFactors = F,sep="\t",header=F)[,1]
subl<-sharpats2
pathways<-list()
for(a in 1:length(subl)){
  if(length(subl[[a]]$Downup)>1){
    pat<-gprofiler(subl[[a]]$Downup,organism="hsapiens",custom_bg=background)  ;  pathways[[names(subl)[a]]]$Downup<-pat
    print(paste(names(subl)[a],length(pat$term.name),"down-up",sep="  "))
  }
  if(length(subl[[a]]$Updown)>1){
    pat<-gprofiler(subl[[a]]$Updown,organism="hsapiens",custom_bg=background)  ;  pathways[[names(subl)[a]]]$Updown<-pat
    print(paste(names(subl)[a],length(pat$term.name),"up-down",sep="  "))
  }
}
save(pathways,file="Data/Subgroups_pathways.Rdata")

## Run enrichment
load("Data/Overlapping_pathways.Rdata") #pathways

todas<-list()  ;  genimplicado<-list()
for(a in 1:length(pathways)){
  for(b in 1:length(pathways[[a]])){
    if(length(grep("^NSCLC",names(pathways)[a]))==1){
      if(sum(dim(pathways[[a]][[b]]))>14){
        print(names(pathways)[a])  ;  print(table(pathways[[a]][[b]][,10]))  ;  pati<-pathways[[a]][[b]][,c(12,14)]
        if(sum(dim(pati))>2){
          son<-c()
          for(c in 1:length(pati[,1])){cojo<-which(names(sharpats)==names(pathways)[a])
            if(length(intersect(strsplit(pati[c,2],",")[[1]],sharpats[[cojo]][[which(names(sharpats[[cojo]])==names(pathways[[a]])[b])]]))>0){son<-c(son,pati[c,1])}
          }
          if(names(pathways[[a]][b])=="Updown"){todas[[names(pathways)[a]]]$Updown<-pati[,1]}
          if(names(pathways[[a]][b])=="Downup"){todas[[names(pathways)[a]]]$Downup<-pati[,1]}
          if(names(pathways[[a]][b])=="Updown" && length(son)>0){updown<-son  ;  genimplicado[[names(pathways)[a]]]$Updown<-updown}
          if(names(pathways[[a]][b])=="Downup" && length(son)>0){downup<-son  ;  genimplicado[[names(pathways)[a]]]$Downup<-downup}
        }
      }
    }
  }
}
for(a in 1:length(genimplicado)){
  qui<-strsplit(names(genimplicado)[a],"_")[[1]]
  for(b in 1:length(names(genimplicado[[a]]))){
    if(names(genimplicado[[a]])[b]=="Updown"){proportions[[qui[1]]]$Up  ;  proportions[[qui[2]]]$Down}
    if(names(genimplicado[[a]])[b]=="Downup"){proportions[[qui[1]]]$Down  ;  proportions[[qui[2]]]$Up}
  }
}

ups<-c();downs<-c()
for(a in grep("Alzheimer",names(degs))){
  if(length(which(degs[[a]]$up=="NUP54"))==1){ups<-c(ups,a)}
  if(length(which(degs[[a]]$down=="NUP54"))==1){downs<-c(downs,a)}
}
subs<-info$Cluster_info_cluster_specific  ;  names(subs)<-info$Patient
table(as.character(subs[names(degs)[downs]]))

## In the case for the direct comorbidity between AD and NSCLC subgroups
tabluco<-pathways[[1]][[1]][,c(3,12)]
tabluco[order(tabluco[,1]),]

## Lets generate a network of genes and subgroups 
tabla<-c()
for(a in 1:length(names(sharpats))){
  tienen<-names(sharpats[[a]])
  splited<-strsplit(names(sharpats)[a],"_")[[1]]
  for(b in 1:length(tienen)){
    if(length(sharpats[[a]][[b]])>0){
      if(tienen[b]=="Downup"){tabla<-rbind(tabla,cbind(splited[1],splited[2],sharpats[[a]][[b]],"Down","Up"))}
      if(tienen[b]=="Updown"){tabla<-rbind(tabla,cbind(splited[1],splited[2],sharpats[[a]][[b]],"Up","Down"))}
    }
  }
}
tabla2<-tabla[grep("NSCLC",tabla[,1]),]
net<-c();for(a in 1:length(tabla2[,1])){net<-rbind(net,rbind(c(tabla2[a,1],tabla2[a,2],"-"),c(tabla2[a,1],tabla2[a,3],tabla2[a,4]),c(tabla2[a,2],tabla2[a,3],tabla2[a,5])))}
colnames(net)<-c("Interactor_1","Interactor_2","Interaction")
write.table(net,"Networks/Subgroups_shared_genes_interactions.txt",row.names=F,sep="\t",quote=F)

                                                               ## @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ ##
                                                               ## @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ ##
                                                               ## @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ ##


## Generate the subgroups' network based on the shared pathways ##
## @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @ @ @@ @@ @@ @@ @@ @@ @@ @@ @@ ##
##First run the previous scripts for generating heatmaps
clusterl<-enri_reac
negativast<-cunineg  ;  positivast<-cunipos
negativost<-c()
for(a in 1:length(negativast[,1])){
  unodir<-clusterl[[negativast[a,1]]]$Up  ;  unoinv<-clusterl[[negativast[a,1]]]$Down
  dosdir<-clusterl[[negativast[a,2]]]$Up  ;  dosinv<-clusterl[[negativast[a,2]]]$Down
  dd<-intersect(unodir,dosdir)  ;  ii<-intersect(unoinv,dosinv)  ;  di<-intersect(unodir,dosinv)  ;  id<-intersect(unoinv,dosdir)
  negativost<-rbind(negativost,c(negativast[a,],length(unodir),length(unoinv),length(dosdir),length(dosinv),length(dd),length(ii),length(di),length(id)))
}
colnames(negativost)<-c("Cluster_1","Cluster_2","RR","Direct_1","Inverse_1","Direct_2","Inverse_2","Dir-Dir","Inv-Inv","Dir-Inv","Inv-Dir")
cumplenunoneg<-intersect(which(as.numeric(negativost[,10])>0),intersect(which(as.numeric(negativost[,8])==0),which(as.numeric(negativost[,9])==0)))
cumplendosneg<-intersect(which(as.numeric(negativost[,11])>0),intersect(which(as.numeric(negativost[,8])==0),which(as.numeric(negativost[,9])==0)))
cumplenneg<-unique(c(cumplenunoneg,cumplendosneg))
intnegt<-negativost[cumplenneg,]
dn1<-intersect(which(negativost[,5]=="0"),which(negativost[,4]=="0"))
dn2<-intersect(which(negativost[,7]=="0"),which(negativost[,6]=="0"))
dsneg<-unique(c(dn1,dn2))
positivost<-c()
for(a in 1:length(positivast[,1])){
  unodir<-clusterl[[positivast[a,1]]]$Up  ;  unoinv<-clusterl[[positivast[a,1]]]$Down
  dosdir<-clusterl[[positivast[a,2]]]$Up  ;  dosinv<-clusterl[[positivast[a,2]]]$Down
  dd<-intersect(unodir,dosdir)  ;  ii<-intersect(unoinv,dosinv)  ;  di<-intersect(unodir,dosinv)  ;  id<-intersect(unoinv,dosdir)
  positivost<-rbind(positivost,c(positivast[a,],length(unodir),length(unoinv),length(dosdir),length(dosinv),length(dd),length(ii),length(di),length(id)))
}
colnames(positivost)<-c("Cluster_1","Cluster_2","RR","Direct_1","Inverse_1","Direct_2","Inverse_2","Dir-Dir","Inv-Inv","Dir-Inv","Inv-Dir")
cumplenunopos<-intersect(which(as.numeric(positivost[,8])>0),intersect(which(as.numeric(positivost[,10])==0),which(as.numeric(positivost[,11])==0)))
cumplendospos<-intersect(which(as.numeric(positivost[,9])>0),intersect(which(as.numeric(positivost[,10])==0),which(as.numeric(positivost[,11])==0)))
cumplenpos<-unique(c(cumplenunopos,cumplendospos))
intpost<-positivost[cumplenpos,]
dp1<-intersect(which(positivost[,5]=="0"),which(positivost[,4]=="0"))
dp2<-intersect(which(positivost[,7]=="0"),which(positivost[,6]=="0"))
dspos<-unique(c(dp1,dp2))
## Represent again the heatmap but in this case only with the clusters sharing drugs in the correct direction
interestmapt<-matrix(ncol=length(cenf[,1]),nrow=length(cenf[,1]),0)
rownames(interestmapt)<-rownames(cenf)  ;  colnames(interestmapt)<-colnames(cenf)
for(a in 1:length(intnegt[,1])){interestmapt[as.character(intnegt[a,1]),as.character(intnegt[a,2])]<-as.numeric(intnegt[a,3])*(-1)}
for(a in 1:length(intpost[,1])){interestmapt[as.character(intpost[a,1]),as.character(intpost[a,2])]<-as.numeric(intpost[a,3])}
#Generating the network 
rencilla<-rbind(cbind(intnegt[,1:2],-1),cbind(intpost[,1:2],1))
colnames(rencilla)<-c("Interactor_1","Interactor_2","Interaction")
red<-rbind(cbind(intpost[,c(1,2,3)],1),cbind(intnegt[,c(1,2,3)],-1))
colnames(red)<-c("Interactor_1","Interactor_2","Intensity","Direction")
if(length(which(gsub("\\..+","",red[,1])==gsub("\\..+","",red[,2])))>0){red<-red[-which(gsub("\\..+","",red[,1])==gsub("\\..+","",red[,2])),]}


## Calculate intra-disease interaction percentages for each disease (this is going to be a supplementary table) ##
## @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @ @ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ ##
## Intra-disease and intra-subgroup interaction percentages ##
load("Data/Patient_patient_similarity_network.Rdata")
pats<-unique(c(fisher_table[,1],fisher_table[,2]))
diseases<-unique(gsub("_.+","",pats))
info<-read.csv2("Patients_information_table.txt",stringsAsFactors = F,sep="\t")
clus<-info$Cluster_info_cluster_specific  ;  names(clus)<-info$Patient
perdis<-c()  ;  subperdis<-c()
for(a in 1:length(diseases)){
  # a<-1
  intersac<-intersect(grep(paste("^",diseases[a],"_",sep=""),fisher_table[,1]),grep(paste("^",diseases[a],"_",sep=""),fisher_table[,2]))
  tab<-table(as.character(fisher_table[intersac,3]))
  if(length(which(names(tab)==1))==0){uno<-0}
  if(length(which(names(tab)==1))>0){uno<-as.character(tab[which(names(tab)==1)])}
  if(length(which(names(tab)==0))==0){cero<-0}
  if(length(which(names(tab)==0))>0){cero<-as.character(tab[which(names(tab)==0)])}
  if(length(which(names(tab)==-1))==0){dos<-0}
  if(length(which(names(tab)==-1))>0){dos<-as.character(tab[which(names(tab)==-1)])}
  perdis<-rbind(perdis,c(diseases[a],uno,dos,cero))
  peque<-cbind(fisher_table[intersac,1],fisher_table[intersac,2],as.character(fisher_table[intersac,3]))
  uns<-as.character(clus[peque[,1]])  ;  uds<-as.character(clus[peque[,2]])
  subs<-unique(c(uns,uds))
  for(b in 1:length(subs)){
    intersac<-intersect(which(uns==subs[b]),which(uds==subs[b]))
    tab<-table(as.character(peque[intersac,3]))
    if(length(which(names(tab)==1))==0){uno<-0}
    if(length(which(names(tab)==1))>0){uno<-as.character(tab[which(names(tab)==1)])}
    if(length(which(names(tab)==0))==0){cero<-0}
    if(length(which(names(tab)==0))>0){cero<-as.character(tab[which(names(tab)==0)])}
    if(length(which(names(tab)==-1))==0){dos<-0}
    if(length(which(names(tab)==-1))>0){dos<-as.character(tab[which(names(tab)==-1)])}
    subperdis<-rbind(subperdis,c(subs[b],uno,dos,cero))
  }
  print(a)
}
colnames(perdis)<-c("Disease","Positivo","Negativo","Cero")
porcen<-as.numeric(perdis[,2])/(as.numeric(perdis[,2])+as.numeric(perdis[,3])+as.numeric(perdis[,4]))
perdis2<-cbind(perdis,porcen)
perdis2[order(porcen,decreasing=T),]

## Drugs by subgroup, suffled vs. true ##
## @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ ##
## We do this only for the true subgroups (those presenting more shared genes than when assigning randomly patients from the same disease to subgroups)
set.seed(1)
info<-read.csv2("Patients_information_table.txt",stringsAsFactors = F,sep="\t")
patients_subgroups<-info[,c(1,2,8)]
subgrupos<-unique(patients_subgroups[,3])  ;  subgrupos<-subgrupos[-c(grep("^^",subgrupos,fixed=T),grep("NS$",subgrupos))]
todos<-patients_subgroups[,c(1,3)]  ;  todos<-todos[-grep("^^",todos[,2],fixed=T),]
load("Data/Drugs.Rdata") # drugs ; load drugs associated to each patient
print("We start shuffling")
totales<-c()
listorras<-list()
for(z in 1:10000){
  drogacom<-c()  ;  drogascomunes<-list()
  for(a in 1:length(subgrupos)){
    coge<-subgrupos[a]
    ## randomly assign patients from the same disease to subgroups
    tos<-todos[grep(strsplit(coge,".",fixed=T)[[1]][1],todos[,2]),]  ;  todos2<-cbind(tos[,1],sample(tos[,2]))
    pacientes<-todos2[which(todos2[,2]==coge),1]  ;  similarl<-list();oppositel<-list()
    for(b in 1:length(pacientes)){similarl[[pacientes[b]]]<-drugs[pacientes[b]][[1]]$Similar  ;  oppositel[[pacientes[b]]]<-drugs[pacientes[b]][[1]]$Opposite}
    similars<-Reduce(intersect,similarl)  ;  opposites<-Reduce(intersect,oppositel)
    drogascomunes[[subgrupos[a]]]$Similar<-similars  ;  drogascomunes[[subgrupos[a]]]$Opposite<-opposites
    if(length(similars)>0 || length(opposites)>0){drogacom<-c(drogacom,1)}
  }
  totales<-c(totales,length(drogacom))  ;  print(paste(z,length(drogacom),sep="  "))  ;  listorras[[paste("Random_",z,sep="")]]<-drogascomunes
} ## repeat this with a different seed
save(listorras,file="Data/Subgroups_common_drugs_shuffled_intra_disease.Rdata")  ;  print("Finished!!")
save(totales, file="Data/Number_of_drugs_associated_to_shuffled_patients.Rdata")

load("Data/Subgroups_common_drugs_shuffled_intra_disease.Rdata")
load("Data/Number_of_drugs_associated_to_shuffled_patients.Rdata")
## Extract the real number of drugs per subgroup ##
clusterl<-list()  ;  drogacom<-c()
for(a in 1:length(subgrupos)){
  patients<-info$Patient[which(info$Cluster_info_cluster_specific==subgrupos[[a]])]
  clus_sim<-list()  ;  clus_opp<-list()
  for(b in 1:length(patients)){
    clus_sim[[b]]<-drugs[[which(names(drugs)==patients[b])]]$Similar
    clus_opp[[b]]<-drugs[[which(names(drugs)==patients[b])]]$Opposite
  }
  clusterl[[subgrupos[a]]]$Similar<-Reduce(intersect,clus_sim)  ;  similars<-Reduce(intersect,clus_sim)
  clusterl[[subgrupos[a]]]$Opposite<-Reduce(intersect,clus_opp)  ;  opposites<-Reduce(intersect,clus_opp)
  if(length(similars)>0 || length(opposites)>0){drogacom<-c(drogacom,1)}
}
#pval
length(which(totales>=length(drogacom)))/length(totales)
plot(density(totales),xlim=c(min(totales)-10,length(drogacom)+5),main="Patient subgroups with shared drugs associated in\nthe same direction (both positively and negatively)",
     cex.main=0.9)
abline(v=length(drogacom),col="red")
pdf(file="Plots/Random_drug_distribution_against_true_results_for_non_NS_subgroups.pdf")
plot(density(totales),xlim=c(min(totales)-10,length(drogacom)+5),main="Patient subgroups with shared drugs associated in\nthe same direction (both positively and negatively)",
     cex.main=0.9)
abline(v=length(drogacom),col="red")
dev.off()


