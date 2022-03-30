## Enrichment analyses ##
library("gProfileR")
load("Data/Clusters_by_shared_genes.Rdata") # clusterl_genes
load("Data/Biogrid_gene_interactors.Rdata") # genl

intneg2<-clusterl_genes$Intnegt 
## Extract the genes involved in the interactions of interest
repes<-c() ; for(a in 1:length(intneg2[,1])){repes<-rbind(repes,sort(intneg2[a,1:2]))}
intneg<-intneg2[-which(duplicated(repes)),]

sharpats<-list()
for(a in 1:length(intneg[,1])){
  genis<-c(intersect(clusterl_genes$Clusterl[[which(names(clusterl_genes$Clusterl)==intneg[a,1])]]$Up_Dir,
                     clusterl_genes$Clusterl[[which(names(clusterl_genes$Clusterl)==intneg[a,2])]]$Down_Inv),
           intersect(clusterl_genes$Clusterl[[which(names(clusterl_genes$Clusterl)==intneg[a,1])]]$Down_Inv,
                     clusterl_genes$Clusterl[[which(names(clusterl_genes$Clusterl)==intneg[a,2])]]$Up_Dir))
  print(paste(intneg[a,1]," ",intneg[a,2]," -> ",genis,sep=""))
  updown<-intersect(clusterl_genes$Clusterl[[which(names(clusterl_genes$Clusterl)==intneg[a,1])]]$Up_Dir,
                    clusterl_genes$Clusterl[[which(names(clusterl_genes$Clusterl)==intneg[a,2])]]$Down_Inv)
  downup<-intersect(clusterl_genes$Clusterl[[which(names(clusterl_genes$Clusterl)==intneg[a,1])]]$Down_Inv,
                    clusterl_genes$Clusterl[[which(names(clusterl_genes$Clusterl)==intneg[a,2])]]$Up_Dir)
  sharpats[[paste(intneg[a,1],intneg[a,2],sep="_")]]$Updown<-updown
  sharpats[[paste(intneg[a,1],intneg[a,2],sep="_")]]$Downup<-downup
}
sharpatsi<-sharpats

#Generate a list with all the subgroups' genes and their first neighbours for enrichment analyses
sharpats2<-list()
for(a in 1:length(names(sharpats))){
  upi<-c()  ;  downi<-c()
  if(length(sharpats[[a]]$Updown)>0){mis<-sharpats[[a]]$Updown  ;  ups<-c()  ;  for(b in 1:length(mis)){ups<-c(ups,genl[[mis[b]]])}  ;  upi<-unique(c(ups,mis))  ;  sharpats2[[names(sharpats)[a]]]$Updown<-upi}
  if(length(upi)==0){sharpats2[[names(sharpats)[a]]]$Updown<-""}
  if(length(sharpats[[a]]$Downup)>0){mis<-sharpats[[a]]$Downup  ;  downs<-c();for(b in 1:length(mis)){downs<-c(downs,genl[[mis[b]]])}
  downi<-unique(c(downs,mis))  ;  sharpats2[[names(sharpats)[a]]]$Downup<-downi}  ;  if(length(downi)==0){sharpats2[[names(sharpats)[a]]]$Downup<-""}
}

save(sharpats2,file="Data/Enrichments/Shared_genes_inverse_comorbidity_size_4_subgroups.Rdata")

## Run enrichment on inverse comorbidities ##
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
pathwaysi<-pathways
save(pathwaysi,file="Data/Enrichments/Subgroups_pathways_inv.Rdata")

## Positives ##
## @@ @ @ @@ ##
## Enrichment analyses ##
library("gProfileR")
load("Data/Clusters_by_shared_genes.Rdata") # clusterl_genes
load("Data/Biogrid_gene_interactors.Rdata") # genl
intpos<-clusterl_genes$Intpost ; intpos2<-intpos[which(gsub("\\..+","",intpos[,1])!=gsub("\\..+","",intpos[,2])),]
## Extract the genes involved in the interactions of interest
repes<-c() ; for(a in 1:length(intpos2[,1])){repes<-rbind(repes,sort(intpos2[a,1:2]))}
intpos<-intpos2[-which(duplicated(repes)),]
## Extract the genes involved in the interactions of interest
sharpats<-list()
for(a in 1:length(intpos[,1])){
  genis<-c(intersect(clusterl_genes$Clusterl[[which(names(clusterl_genes$Clusterl)==intpos[a,1])]]$Up_Dir,
                     clusterl_genes$Clusterl[[which(names(clusterl_genes$Clusterl)==intpos[a,2])]]$Up_Dir),
           intersect(clusterl_genes$Clusterl[[which(names(clusterl_genes$Clusterl)==intpos[a,1])]]$Down_Inv,
                     clusterl_genes$Clusterl[[which(names(clusterl_genes$Clusterl)==intpos[a,2])]]$Down_Inv))
  print(paste(intpos[a,1]," ",intpos[a,2]," -> ",genis,sep=""))
  updown<-intersect(clusterl_genes$Clusterl[[which(names(clusterl_genes$Clusterl)==intpos[a,1])]]$Up_Dir,
                    clusterl_genes$Clusterl[[which(names(clusterl_genes$Clusterl)==intpos[a,2])]]$Up_Dir)
  downup<-intersect(clusterl_genes$Clusterl[[which(names(clusterl_genes$Clusterl)==intpos[a,1])]]$Down_Inv,
                    clusterl_genes$Clusterl[[which(names(clusterl_genes$Clusterl)==intpos[a,2])]]$Down_Inv)
  sharpats[[paste(intpos[a,1],intpos[a,2],sep="_")]]$UpUp<-updown
  sharpats[[paste(intpos[a,1],intpos[a,2],sep="_")]]$DownDown<-downup
}
sharpatsd<-sharpats
#Generate a list with all the subgroups' genes and their first neighbours for enrichment analyses
sharpats2<-list()
for(a in 1:length(names(sharpats))){
  upi<-c()  ;  downi<-c()
  if(length(sharpats[[a]]$UpUp)>0){mis<-sharpats[[a]]$UpUp  ;  ups<-c()  ;  for(b in 1:length(mis)){ups<-c(ups,genl[[mis[b]]])}  ;  upi<-unique(c(ups,mis))  ;  sharpats2[[names(sharpats)[a]]]$UpUp<-upi}
  if(length(upi)==0){sharpats2[[names(sharpats)[a]]]$UpUp<-""}
  if(length(sharpats[[a]]$DownDown)>0){mis<-sharpats[[a]]$DownDown  ;  downs<-c();for(b in 1:length(mis)){downs<-c(downs,genl[[mis[b]]])}
  downi<-unique(c(downs,mis))  ;  sharpats2[[names(sharpats)[a]]]$DownDown<-downi}  ;  if(length(downi)==0){sharpats2[[names(sharpats)[a]]]$DownDown<-""}
}
save(sharpats2,file="Data/Enrichments/Shared_genes_direct_comorbidity_size_4_subgroups.Rdata")
## Run enrichment on inverse comorbidities ##
background<-read.csv2("Data/Background.txt",stringsAsFactors = F,sep="\t",header=F)[,1]
subl<-sharpats2
pathways<-list()
for(a in 1:length(subl)){
  if(length(subl[[a]]$DownDown)>1){
    pat<-gprofiler(subl[[a]]$DownDown,organism="hsapiens",custom_bg=background)  ;  pathways[[names(subl)[a]]]$DownDown<-pat
    print(paste(names(subl)[a],length(pat$term.name),"down-down",sep="  "))
  }
  if(length(subl[[a]]$UpUp)>1){
    pat<-gprofiler(subl[[a]]$UpUp,organism="hsapiens",custom_bg=background)  ;  pathways[[names(subl)[a]]]$UpUp<-pat
    print(paste(names(subl)[a],length(pat$term.name),"up-up",sep="  "))
  }
}

pathwaysd<-pathways
save(pathwaysd,file="Data/Enrichments/Subgroups_pathways_dir.Rdata")
