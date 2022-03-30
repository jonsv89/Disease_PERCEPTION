setwd("./")
## Extract identifiers and then calculate Relative Risk between ICD9s, ICD10s, diseases & patient-subgroups ##
## Confidence intervals ##
## 95% --> 1.96
## 99% --> 2.576
## load function
estimate_rr<-function(disease,identificadores,inter0001,path,confidenceinterval){
  RRpos<-c()
  RRneg<-c()
  for(a in 1:(length(disease)-1)){
    beginloop<-Sys.time()
    ## calculate relative risks for each disease 
    print(disease[a])
    ## select first diseases' interactions of interest
    enf1<-identificadores[[disease[a]]]
    for(b in (a+1):length(disease)){
      notpresentbeing<-c() ; presentbeing<-c() ; presentnotbeing<-c() ; notpresentnotbeing<-c() ; notpresentbeing1<-c() ; notpresentbeing2<-c()
      notpresentnotbeing1<-c() ; notpresentnotbeing2<-c()
      ## select second diseases' interactions of interest
      enf2<-identificadores[[disease[b]]]
      par<-intersect(enf1,enf2)
      only1<-setdiff(enf1,par)
      only2<-setdiff(enf2,par)
      
      ## @ @ ##
      ## pRR ##
      ## @ @ ##
      being<-table(inter0001[par,3])
      # Not present being
      if(length(which(names(being)==0))>0){notpresentbeing1<-as.numeric(being[which(names(being)==0)])}
      if(length(which(names(being)==(-1)))>0){notpresentbeing2<-as.numeric(being[which(names(being)==(-1))])}
      if(length(notpresentbeing1)>0){notpresentbeing<-notpresentbeing1}
      if(length(notpresentbeing2)>0){notpresentbeing<-notpresentbeing2}
      if(length(notpresentbeing1)>0 && length(notpresentbeing2)>0){notpresentbeing<-notpresentbeing1+notpresentbeing2}
      if(length(notpresentbeing)==0){notpresentbeing<-0}
      # Present being
      if(length(which(names(being)==1))>0){presentbeing<-as.numeric(being[which(names(being)==1)])}
      if(length(presentbeing)==0){presentbeing<-0}
      ## Disease 1 ##
      ## @ @@ @@ @ ##
      notbeing<-table(inter0001[only1,3])
      # Not present not being
      if(length(which(names(notbeing)==0))>0){notpresentnotbeing1<-as.numeric(notbeing[which(names(notbeing)==0)])}
      if(length(which(names(notbeing)==(-1)))>0){notpresentnotbeing2<-as.numeric(notbeing[which(names(notbeing)==(-1))])}
      if(length(notpresentnotbeing1)>0){notpresentnotbeing<-notpresentnotbeing1}
      if(length(notpresentnotbeing2)>0){notpresentnotbeing<-notpresentnotbeing2}
      if(length(notpresentnotbeing1)>0 && length(notpresentnotbeing2)>0){notpresentnotbeing<-notpresentnotbeing1+notpresentnotbeing2}
      if(length(notpresentnotbeing)==0){notpresentnotbeing<-0}
      # Present not being
      if(length(which(names(notbeing)==1))>0){presentnotbeing<-as.numeric(notbeing[which(names(notbeing)==1)])}
      if(length(presentnotbeing)==0){presentnotbeing<-0}
      peventwhenexposed<-presentbeing/(presentbeing+notpresentbeing)
      peventwhennotexposed<-presentnotbeing/(presentnotbeing+notpresentnotbeing)
      division<-peventwhenexposed/peventwhennotexposed
      if(length(division)>0){
        raiz<-sqrt(((notpresentbeing/presentbeing)/(notpresentbeing+presentbeing))+((notpresentnotbeing/presentnotbeing)/(notpresentnotbeing+presentnotbeing)))
        sup<-exp(log(division)+(confidenceinterval*raiz))
        inf<-exp(log(division)-(confidenceinterval*raiz))
        RRpos<-rbind(RRpos,c(disease[a],disease[b],division,sup,inf))
        # print(paste(disease[a],disease[b],division,sup,inf,sep="  "))
      }
      ## Disease 2 ##
      ## @ @@ @@ @ ##
      notpresentnotbeing1<-c() ; notpresentnotbeing2<-c() ; presentnotbeing<-c() ; notbeing<-table(inter0001[only2,3])
      # Not present not being
      if(length(which(names(notbeing)==0))>0){notpresentnotbeing1<-as.numeric(notbeing[which(names(notbeing)==0)])}
      if(length(which(names(notbeing)==(-1)))>0){notpresentnotbeing2<-as.numeric(notbeing[which(names(notbeing)==(-1))])}
      if(length(notpresentnotbeing1)>0){notpresentnotbeing<-notpresentnotbeing1}
      if(length(notpresentnotbeing2)>0){notpresentnotbeing<-notpresentnotbeing2}
      if(length(notpresentnotbeing1)>0 && length(notpresentnotbeing2)>0){notpresentnotbeing<-notpresentnotbeing1+notpresentnotbeing2}
      if(length(notpresentnotbeing)==0){notpresentnotbeing<-0}
      # Present not being
      if(length(which(names(notbeing)==1))>0){presentnotbeing<-as.numeric(notbeing[which(names(notbeing)==1)])}
      if(length(presentnotbeing)==0){presentnotbeing<-0}
      peventwhenexposed<-presentbeing/(presentbeing+notpresentbeing)
      peventwhennotexposed<-presentnotbeing/(presentnotbeing+notpresentnotbeing)
      division<-peventwhenexposed/peventwhennotexposed
      if(length(division)>0){
        raiz<-sqrt(((notpresentbeing/presentbeing)/(notpresentbeing+presentbeing))+((notpresentnotbeing/presentnotbeing)/(notpresentnotbeing+presentnotbeing)))
        sup<-exp(log(division)+(confidenceinterval*raiz))
        inf<-exp(log(division)-(confidenceinterval*raiz))
        RRpos<-rbind(RRpos,c(disease[b],disease[a],division,sup,inf))
        # print(paste(disease[b],disease[a],division,sup,inf,sep="  "))
      }
      
      ## @ @ ##
      ## nRR ##
      ## @ @ ##
      notpresentbeing<-c() ; presentbeing<-c() ; presentnotbeing<-c() ; notpresentnotbeing<-c() ; notpresentbeing1<-c() ; notpresentbeing2<-c()
      notpresentnotbeing1<-c() ; notpresentnotbeing2<-c()
      being<-table(inter0001[par,3])
      # Not present being
      if(length(which(names(being)==0))>0){notpresentbeing1<-as.numeric(being[which(names(being)==0)])}
      if(length(which(names(being)==1))>0){notpresentbeing2<-as.numeric(being[which(names(being)==1)])}
      if(length(notpresentbeing1)>0){notpresentbeing<-notpresentbeing1}
      if(length(notpresentbeing2)>0){notpresentbeing<-notpresentbeing2}
      if(length(notpresentbeing1)>0 && length(notpresentbeing2)>0){notpresentbeing<-notpresentbeing1+notpresentbeing2}
      if(length(notpresentbeing)==0){notpresentbeing<-0}
      # Present being
      if(length(which(names(being)==(-1)))>0){presentbeing<-as.numeric(being[which(names(being)==(-1))])}
      if(length(presentbeing)==0){presentbeing<-0}
      ## Disease 1 ##
      ## @ @@ @@ @ ##
      notbeing<-table(inter0001[only1,3])
      # Not present not being
      if(length(which(names(notbeing)==0))>0){notpresentnotbeing1<-as.numeric(notbeing[which(names(notbeing)==0)])}
      if(length(which(names(notbeing)==1))>0){notpresentnotbeing2<-as.numeric(notbeing[which(names(notbeing)==1)])}
      if(length(notpresentnotbeing1)>0){notpresentnotbeing<-notpresentnotbeing1}
      if(length(notpresentnotbeing2)>0){notpresentnotbeing<-notpresentnotbeing2}
      if(length(notpresentnotbeing1)>0 && length(notpresentnotbeing2)>0){notpresentnotbeing<-notpresentnotbeing1+notpresentnotbeing2}
      if(length(notpresentnotbeing)==0){notpresentnotbeing<-0}
      # Present not being
      if(length(which(names(notbeing)==(-1)))>0){presentnotbeing<-as.numeric(notbeing[which(names(notbeing)==(-1))])}
      if(length(presentnotbeing)==0){presentnotbeing<-0}
      peventwhenexposed<-presentbeing/(presentbeing+notpresentbeing)
      peventwhennotexposed<-presentnotbeing/(presentnotbeing+notpresentnotbeing)
      division<-peventwhenexposed/peventwhennotexposed
      if(length(division)>0){
        raiz<-sqrt(((notpresentbeing/presentbeing)/(notpresentbeing+presentbeing))+((notpresentnotbeing/presentnotbeing)/(notpresentnotbeing+presentnotbeing)))
        sup<-exp(log(division)+(confidenceinterval*raiz))
        inf<-exp(log(division)-(confidenceinterval*raiz))
        RRneg<-rbind(RRneg,c(disease[a],disease[b],division,sup,inf))
        # print(paste(disease[a],disease[b],division,sup,inf,sep="  "))
      }
      ## Disease 2 ##
      ## @ @@ @@ @ ##
      notpresentnotbeing1<-c() ; notpresentnotbeing2<-c() ; presentnotbeing<-c()
      notbeing<-table(inter0001[only2,3])
      # Not present not being
      if(length(which(names(notbeing)==0))>0){notpresentnotbeing1<-as.numeric(notbeing[which(names(notbeing)==0)])}
      if(length(which(names(notbeing)==1))>0){notpresentnotbeing2<-as.numeric(notbeing[which(names(notbeing)==1)])}
      if(length(notpresentnotbeing1)>0){notpresentnotbeing<-notpresentnotbeing1}
      if(length(notpresentnotbeing2)>0){notpresentnotbeing<-notpresentnotbeing2}
      if(length(notpresentnotbeing1)>0 && length(notpresentnotbeing2)>0){notpresentnotbeing<-notpresentnotbeing1+notpresentnotbeing2}
      if(length(notpresentnotbeing)==0){notpresentnotbeing<-0}
      # Present not being
      if(length(which(names(notbeing)==(-1)))>0){presentnotbeing<-as.numeric(notbeing[which(names(notbeing)==(-1))])}
      if(length(presentnotbeing)==0){presentnotbeing<-0}
      peventwhenexposed<-presentbeing/(presentbeing+notpresentbeing)
      peventwhennotexposed<-presentnotbeing/(presentnotbeing+notpresentnotbeing)
      division<-peventwhenexposed/peventwhennotexposed
      if(length(division)>0){
        raiz<-sqrt(((notpresentbeing/presentbeing)/(notpresentbeing+presentbeing))+((notpresentnotbeing/presentnotbeing)/(notpresentnotbeing+presentnotbeing)))
        sup<-exp(log(division)+(confidenceinterval*raiz))
        inf<-exp(log(division)-(confidenceinterval*raiz))
        RRneg<-rbind(RRneg,c(disease[b],disease[a],division,sup,inf))
        # print(paste(disease[b],disease[a],division,sup,inf,sep="  "))
      }
    }
    endloop<-Sys.time()
    save(RRpos,file=paste(path,"pRR.Rdata",sep="_"))
    save(RRneg,file=paste(path,"nRR.Rdata",sep="_"))
    print(paste("Finished ",a," out of ",length(disease)," in ",endloop-beginloop,sep=" "))
  }
  rrs<-list("pRR"=RRpos,"nRR"=RRneg)
  print("Finished!!")
  return(rrs)
}

## generate folders
if("RR"%in%list.files("Data/")==FALSE){dir.create("./Data/RR/")}
if ("Identifiers"%in%list.files("Data/") == FALSE){dir.create("Data/Identifiers")}



## @@ @@ @@ ##
## Diseases ##
## @@ @@ @@ ##
## Generate identifiers ##
## @@ @@ @@ @@ @@ @@ @@ ##
## estimated duration: 34.12 mins
beggining<-Sys.time()
load("Data/Networks_varying_thresholds/Patient_similarities/Patient_patient_similarity_network_1e-04.Rdata") # fisher_table
# load("Data/Patient_patient_similarity_network.Rdata") # fisher_table
net<-fisher_table  ;  net0001<-net  ;  disease<-unique(gsub("_.+","",unique(c(net0001[,1],net0001[,2]))))
one<-gsub("_.+","",net0001[,1],fixed=F)  ;  two<-gsub("_.+","",net0001[,2],fixed=F)
intra<-which(one==two)  ;  inter0001<-net0001[-intra,]  ;  identificadores<-list()
print("Starting the loop")
for(a in 1:length(disease)){
  ini<-Sys.time();identificadores[[disease[a]]]<-c(grep(paste(disease[a],"_",sep=""),inter0001[,1]),grep(paste(disease[a],"_",sep=""),inter0001[,2]));fin<-Sys.time()
  print(paste(a," out of ",length(disease)," in ",fin-ini,sep=""))
}
save(identificadores,file="Data/Identifiers/Disease.Rdata")
ending<-Sys.time()
print(paste("Finished in: ",ending-beggining,"!",sep=""))
## Calculate Relative Risks ##
## @@ @@ @@ @@  @@ @@ @@ @@ ##
## estimated duration: 12.35 minutes
starting<-Sys.time()
load("Data/Networks_varying_thresholds/Patient_similarities/Patient_patient_similarity_network_1e-04.Rdata")
net<-fisher_table
net0001<-net
disease<-unique(gsub("_.+","",unique(c(net0001[,1],net0001[,2]))))
## remove intra-disease interactions
one<-gsub("_.+","",net0001[,1],fixed=F)  ;  two<-gsub("_.+","",net0001[,2],fixed=F)  ;  intra<-which(one==two)  ;  inter0001<-net0001[-intra,]
## load previously generated identifiers
load("Data/Identifiers/Disease.Rdata")
path<-"./Data/RR/Disease"
estimate_rr(disease,identificadores,inter0001,path,1.96)
finishing<-Sys.time()  ;  print(finishing-starting)



## @@ @@ @@
## ICD9s ##
## @@ @@ ##
## Generate identifiers ##
## @@ @@ @@ @@ @@ @@ @@ ##
## estimated duration: 20.06 mins
beggining<-Sys.time()
load("Data/Networks_varying_thresholds/Patient_similarities/Patient_patient_similarity_network_1e-04.Rdata") # fisher_table
# load("Data/Patient_patient_similarity_network.Rdata") # fisher_table
## do the same as in the last step, just changing disease name by the corresponding ICD9 code
patinf<-read.csv2("Patients_information_table.txt",stringsAsFactors = F,sep="\t",header=T)
splited<-strsplit(patinf[,1],"_")  ;  newnames<-c()
for(a in 1:length(splited)){newnames<-c(newnames,paste(gsub("_","-",patinf[a,4]),splited[[a]][2],splited[[a]][3],splited[[a]][4],sep="_"))}
names(newnames)<-patinf[,1]
net<-cbind(as.character(newnames[as.character(fisher_table[,1])]),as.character(newnames[as.character(fisher_table[,2])]),as.character(fisher_table[,3]))
colnames(net)<-colnames(fisher_table)  ;  net0001<-net  ;  disease<-unique(gsub("_.+","",unique(c(net0001[,1],net0001[,2]))))
one<-gsub("_.+","",net0001[,1],fixed=F)  ;  two<-gsub("_.+","",net0001[,2],fixed=F)  ;  intra<-which(one==two)  ;  inter0001<-net0001[-intra,]
identificadores<-list()  ;  print("Starting the loop")
for(a in 1:length(disease)){
  ini<-Sys.time();identificadores[[disease[a]]]<-c(grep(paste(disease[a],"_",sep=""),inter0001[,1]),grep(paste(disease[a],"_",sep=""),inter0001[,2]));fin<-Sys.time()
  print(paste(a," out of ",length(disease)," in ",fin-ini,sep=""))
}
save(identificadores,file="Data/Identifiers/ICD9.Rdata")
ending<-Sys.time()  ;  print(paste("Finished in: ",ending-beggining,"!",sep=""))
## Calculate Relative Risks ##
## @@ @@ @@ @@  @@ @@ @@ @@ ##
## estimated duration: 10.68 minutes
starting<-Sys.time()
## load previously generated identifiers
load("Data/Identifiers/ICD9.Rdata")
path<-"./Data/RR/ICD9"
## since this network is being generated to compare it with the PDN we use a more restrictive threshold (the one used by Hidalgo et al. 2009)
estimate_rr(disease,identificadores,inter0001,path,2.576)
finishing<-Sys.time()  ;  print(finishing-starting)



## @@  @@ ##
## ICD10s ##
## @@  @@ ##
## Generate identifiers ##
## @@ @@ @@ @@ @@ @@ @@ ##
## estimated duration: 20.06 mins
beggining<-Sys.time()
load("Data/Networks_varying_thresholds/Patient_similarities/Patient_patient_similarity_network_1e-04.Rdata") # fisher_table
# load("Data/Patient_patient_similarity_network.Rdata") # fisher_table
## do the same as in the last step, just changing disease name by the corresponding ICD10 code
patinf<-read.csv2("Patients_information_table.txt",stringsAsFactors = F,sep="\t",header=T)
splited<-strsplit(patinf[,1],"_")  ;  newnames<-c()
for(a in 1:length(splited)){newnames<-c(newnames,paste(patinf[a,11],splited[[a]][2],splited[[a]][3],splited[[a]][4],sep="_"))}
names(newnames)<-patinf[,1]
net<-cbind(as.character(newnames[as.character(fisher_table[,1])]),as.character(newnames[as.character(fisher_table[,2])]),as.character(fisher_table[,3]))
colnames(net)<-colnames(fisher_table)  ;  net0001<-net  ;  disease<-unique(gsub("_.+","",unique(c(net0001[,1],net0001[,2]))))
one<-gsub("_.+","",net0001[,1],fixed=F)  ;  two<-gsub("_.+","",net0001[,2],fixed=F)  ;  intra<-which(one==two)  ;  inter0001<-net0001[-intra,]
identificadores<-list()  ;  print("Starting the loop")
for(a in 1:length(disease)){
  ini<-Sys.time();identificadores[[disease[a]]]<-c(grep(paste(disease[a],"_",sep=""),inter0001[,1]),grep(paste(disease[a],"_",sep=""),inter0001[,2]));fin<-Sys.time()
  print(paste(a," out of ",length(disease)," in ",fin-ini,sep=""))
}
save(identificadores,file="Data/Identifiers/ICD10.Rdata")
ending<-Sys.time()  ;  print(paste("Finished in: ",ending-beggining,"!",sep=""))
## Calculate Relative Risks ##
## @@ @@ @@ @@  @@ @@ @@ @@ ##
## estimated duration: 11.54 minutes
starting<-Sys.time()
## remove intra-disease interactions
load("Data/Identifiers/ICD10.Rdata")
path<-"./Data/RR/ICD10"
estimate_rr(disease,identificadores,inter0001,path,1.96)
finishing<-Sys.time()  ;  print(finishing-starting)



## @@ @@ @@ @@ @@ @@ ##
## Patient-subgroups ##
## @@ @@ @@ @@ @@ @@ ##
## Generate identifiers ##
## @@ @@ @@ @@ @@ @@ @@ ##
## estimated duration: 5.54 hours
beggining<-Sys.time()
load("Data/Networks_varying_thresholds/Patient_similarities/Patient_patient_similarity_network_1e-04.Rdata") # fisher_table
# load("Data/Patient_patient_similarity_network.Rdata") # fisher_table
## do the same as in the last step, just changing disease name by the corresponding ICD9 code
patinf<-read.csv2("Patients_information_table.txt",stringsAsFactors = F,sep="\t",header=T)
splited<-strsplit(patinf[,1],"_")  ;  newnames<-c()
for(a in 1:length(splited)){newnames<-c(newnames,paste(gsub("_","-",patinf[a,8]),splited[[a]][2],splited[[a]][3],splited[[a]][4],sep="_"))}  ;  names(newnames)<-patinf[,1]
net<-cbind(as.character(newnames[as.character(fisher_table[,1])]),as.character(newnames[as.character(fisher_table[,2])]),as.character(fisher_table[,3]))
colnames(net)<-colnames(fisher_table)  ;  net0001<-net
disease<-unique(gsub("_.+","",unique(c(net0001[,1],net0001[,2]))))
one<-gsub("_.+","",net0001[,1],fixed=F)  ;  two<-gsub("_.+","",net0001[,2],fixed=F)  ;  intra<-which(one==two)  ;  inter0001<-net0001[-intra,]
identificadores<-list()  ;  print("Starting the loop")
for(a in 1:length(disease)){
  ini<-Sys.time();identificadores[[disease[a]]]<-c(grep(paste(disease[a],"_",sep=""),inter0001[,1]),grep(paste(disease[a],"_",sep=""),inter0001[,2]));fin<-Sys.time()
  print(paste(a," out of ",length(disease)," in ",fin-ini,sep=""))
}
save(identificadores,file="Data/Identifiers/Patient_subgroup.Rdata")
ending<-Sys.time()  ;  print(paste("Finished in: ",ending-beggining,"!",sep=""))
## Calculate Relative Risks ##
## @@ @@ @@ @@  @@ @@ @@ @@ ##
## estimated duration: 12.35 minutes
starting<-Sys.time()
load("Data/Networks_varying_thresholds/Patient_similarities/Patient_patient_similarity_network_1e-04.Rdata") # fisher_table
patinf<-read.csv2("Patients_information_table.txt",stringsAsFactors = F,sep="\t",header=T)  ;  splited<-strsplit(patinf[,1],"_")  ;  newnames<-c()
for(a in 1:length(splited)){newnames<-c(newnames,paste(gsub("_","-",patinf[a,8]),splited[[a]][2],splited[[a]][3],splited[[a]][4],sep="_"))}
names(newnames)<-patinf[,1]
net<-cbind(as.character(newnames[as.character(fisher_table[,1])]),as.character(newnames[as.character(fisher_table[,2])]),as.character(fisher_table[,3]))
colnames(net)<-colnames(fisher_table)  ;  net0001<-net  ;  disease<-unique(gsub("_.+","",unique(c(net0001[,1],net0001[,2]))))
## remove intra-disease interactions
one<-gsub("_.+","",net0001[,1],fixed=F)  ;  two<-gsub("_.+","",net0001[,2],fixed=F)  ;  intra<-which(one==two)  ;  inter0001<-net0001[-intra,]
## load previously generated identifiers
load("Data/Identifiers/Patient_subgroup.Rdata")
path<-"./Data/RR/Patient_subgroup"
estimate_rr(disease,identificadores,inter0001,path,1.96)
finishing<-Sys.time()  ;  print(finishing-starting)


## @@ @@ @@ @@  @@ @@ @@ @@ ##
## Evaluate network overlap ##
## @@ @@ @@ @@  @@ @@ @@ @@ ##
if("Network_analysis"%in%list.files("Data/")==FALSE){dir.create("Data/Network_analysis")}

set.seed(1)
# dtp<-network_overlap(RRpos,RRneg,network,"Plots/Pairs_DTP-Molecular_Disease_Similarity_Network.pdf")

## load function
network_overlap<-function(RRpos,RRneg,network,output){
  # Select only those codes analyzed in both diseases
  uno<-as.character(network[,1]);dos<-as.character(network[,2])
  interactores<-intersect(unique(c(RRpos[,1],RRpos[,2],RRneg[,1],RRneg[,2])),unique(c(uno,dos)))
  indu<-c()   ;   indd<-c();for(a in 1:length(interactores)){indu<-c(indu,which(uno==interactores[a]))   ;   indd<-c(indd,which(dos==interactores[a]))}
  network2<-cbind(uno[intersect(indu,indd)],dos[intersect(indu,indd)],network[intersect(indu,indd),3])
  network2<-network2[which(as.numeric(network2[,3])>1),]
  print(paste((length(unique(c(network2[,1],network2[,2])))/length(interactores))*100,"% of the diseases interact between them",sep=""))
  notconnected<-setdiff(interactores,unique(c(network2[,1],network2[,2])))
  print(paste(length(unique(c(network2[,1],network2[,2])))," out of ",length(interactores),sep=""))
  # Select significant positive and negative RR interactions, removing those in common
  posit<-RRpos[which(as.numeric(RRpos[,5])>1),1:2]   ;   negat<-RRneg[which(as.numeric(RRneg[,5])>1),1:2]
  indu<-c()   ;   indd<-c();for(a in 1:length(interactores)){indu<-c(indu,which(posit[,1]==interactores[a]))   ;   indd<-c(indd,which(posit[,2]==interactores[a]))}
  positivos<-posit[intersect(indu,indd),]
  indu<-c()   ;   indd<-c();for(a in 1:length(interactores)){indu<-c(indu,which(negat[,1]==interactores[a]))   ;   indd<-c(indd,which(negat[,2]==interactores[a]))}
  negativos<-negat[intersect(indu,indd),]
  posipeg<-paste(positivos[,1],positivos[,2],sep="_")   ;   negapeg<-paste(negativos[,1],negativos[,2],sep="_")
  posipeg2<-setdiff(posipeg,negapeg)   ;   negapeg2<-setdiff(negapeg,posipeg)   ;   posipeg<-posipeg2   ;   negapeg<-negapeg2
  netwpeg<-paste(network2[,1],network2[,2],sep="_")
  # remove next row
  totalnumbers<-c()
  if(length(notconnected)>0){
    totalnumbers<-c();for(buf in 1:length(notconnected)){totalnumbers<-rbind(totalnumbers,c(notconnected[buf],length(grep(notconnected[buf],posipeg))))}
    totalnumbers<-totalnumbers[order(as.numeric(totalnumbers[,2]),decreasing=T),]
  }

  # Generate all the possible interactions between the commonly analyzed ICD codes 
  nueva<-c();for(a in 1:length(interactores)){for(b in 1:length(interactores)){nueva<-rbind(nueva,c(interactores[a],interactores[b]))}}
  nueva<-nueva[-which(nueva[,1]==nueva[,2]),];pasnu<-paste(nueva[,1],nueva[,2],sep="_")
  # Shuffle interactions and calculate overlaps between our positive and negative interactions with the shuffled ones
  ini<-Sys.time()   ;   cuneg<-c()   ;   cunpos<-c()
  for(a in 1:100000){
    aleatorio<-pasnu[sample(1:length(pasnu),length(posipeg))]
    cunpos<-c(cunpos,length(intersect(aleatorio,netwpeg)))
    cuneg<-c(cuneg,length(intersect(aleatorio,netwpeg)))
  }
  print(Sys.time()-ini)
  # Significance of the overlap
  print(length(which(cunpos>=length(intersect(posipeg,netwpeg))))/length(cunpos))
  print(length(which(cuneg>=length(intersect(negapeg,netwpeg))))/length(cuneg))
  pdf(file=output)
  if(length(intersect(posipeg,netwpeg))>=max(cunpos)){
    plot(density(cunpos), main="PDN - pRR interaction intersection",xlim=c(min(cunpos),length(intersect(posipeg,netwpeg))+5))
    abline(v=length(intersect(posipeg,netwpeg)),col="red")
  }
  if(length(intersect(posipeg,netwpeg))<=min(cunpos)){
    plot(density(cunpos), main="PDN - pRR interaction intersection",xlim=c(length(intersect(posipeg,netwpeg))-5,max(cunpos)))
    abline(v=length(intersect(posipeg,netwpeg)),col="red")
  }
  if(length(intersect(posipeg,netwpeg))>=min(cunpos) && length(intersect(posipeg,netwpeg))<=max(cunpos)){
    plot(density(cunpos), main="PDN - pRR interaction intersection")
    abline(v=length(intersect(posipeg,netwpeg)),col="red")
    hist(cunpos,breaks=unique(cunpos),main="PDN - pRR interaction intersection")
    abline(v=length(intersect(posipeg,netwpeg)),col="red")
  }
  if(length(intersect(negapeg,netwpeg))>=max(cuneg)){
    plot(density(cuneg), main="PDN - nRR interaction intersection",xlim=c(min(cuneg),length(intersect(negapeg,netwpeg))+5))
    abline(v=length(intersect(negapeg,netwpeg)),col="red")
  }
  if(length(intersect(negapeg,netwpeg))<=min(cuneg)){
    plot(density(cuneg), main="PDN - nRR interaction intersection",xlim=c(length(intersect(negapeg,netwpeg))-5,max(cuneg)))
    abline(v=length(intersect(negapeg,netwpeg)),col="red")
  }
  if(length(intersect(negapeg,netwpeg))>=min(cuneg) && length(intersect(negapeg,netwpeg))<=max(cuneg)){
    plot(density(cuneg), main="PDN - nRR interaction intersection")
    abline(v=length(intersect(negapeg,netwpeg)),col="red")
    hist(cunpos,breaks=unique(cuneg),main="PDN - nRR interaction intersection")
    abline(v=length(intersect(negapeg,netwpeg)),col="red")
  }
  dev.off()
  results<-list("Epidemiological"=netwpeg,"Positive"=posipeg,"Negative"=negapeg,"Notconnected"=notconnected,"TotalNumbers"=totalnumbers)
  return(results)
}
info<-read.csv2("Patients_information_table.txt",stringsAsFactors = F,sep="\t",header=T)
asoc<-info[,c(3,4,11)][-which(duplicated(info[,c(4,11)])),]

## Disease Trajectories - DMSN ##
## @@ @@ @@ @@ @ @ @@ @@ @@ @@ ##
# Load our ICD10 interaction network and the disease trajectories, selecting only the significant interactions
load("Data/RR/ICD10_pRR.Rdata")
load("Data/RR/ICD10_nRR.Rdata")
network<-read.csv2("Epidemiological_networks/Disease_pairs.csv",stringsAsFactors = F,sep=",",header=T)[,1:3]
dtp<-network_overlap(RRpos,RRneg,network,"Plots/Pairs_DTP-Molecular_Disease_Similarity_Network.pdf")

cates<-info$Disease_group_name
names(cates)<-info$ICD10
tt<-cbind(as.character(cates[gsub("_.+","",dtp$Epidemiological)]),as.character(cates[gsub(".+_","",dtp$Epidemiological)]))
epiicd10int<-as.data.frame(table(paste(tt[,1],tt[,2],sep=" --> ")))[order(as.data.frame(table(paste(tt[,1],tt[,2],sep=" --> ")))[,2],decreasing=T),]
write.table(epiicd10int,"Results/ICD10_comorbidities_categories.txt",quote=F,sep="\t",row.names=F,col.names=F)
## Save the different networks
complete<-network[,1:2]
write.table(complete,"Data/Network_analysis/Disease_trajectories_complete.txt",quote=F,sep="\t",col.names=F,row.names=F)
dtpepi<-cbind(gsub("_.+","",dtp$Epidemiological),gsub(".+_","",dtp$Epidemiological))
write.table(dtpepi,"Data/Network_analysis/Disease_trajectories_subset.txt",quote=F,sep="\t",col.names=F,row.names=F)
dtppos<-cbind(gsub("_.+","",dtp$Positive),gsub(".+_","",dtp$Positive))
write.table(dtppos,"Data/Network_analysis/ICD10_positive.txt",quote=F,sep="\t",col.names=F,row.names=F)
dtpneg<-cbind(gsub("_.+","",dtp$Negative),gsub(".+_","",dtp$Negative))
write.table(dtpneg,"Data/Network_analysis/ICD10_negative.txt",quote=F,sep="\t",col.names=F,row.names=F)
dtp<-dtp$Epidemiological
save(dtp,file="Data/DTP_shared_diseases.Rdata")

## @@ @@@@ @@ ##
## PDN - DMSN ##
## both directions
## @@ @@ @ @ @@ @@
# Load our ICD9 interaction network and the PDN, selecting only the significant interactions
load("Data/RR/ICD9_pRR.Rdata")
load("Data/RR/ICD9_nRR.Rdata")
husig<-read.csv2("Data/PDN_hidalgo.net",stringsAsFactors = F,sep="\t",header=F)   
network2<-cbind(paste("ICD9-",husig[,1],sep=""),paste("ICD9-",husig[,2],sep=""),husig[,7])
network<-rbind(cbind(as.character(network2[,1]),as.character(network2[,2]),network2[,3]),
               cbind(as.character(network2[,2]),as.character(network2[,1]),network2[,3]))
pdn2<-network_overlap(RRpos,RRneg,network,"Plots/PDN-Molecular_Disease_Similarity_Network_both_directions.pdf")

cates<-info$Disease_group_name
names(cates)<-gsub("_","-",info$ICD9)
tt<-cbind(as.character(cates[gsub("_.+","",pdn2$Epidemiological)]),as.character(cates[gsub(".+_","",pdn2$Epidemiological)]))
epiicd9int<-as.data.frame(table(paste(tt[,1],tt[,2],sep=" <--> ")))[order(as.data.frame(table(paste(tt[,1],tt[,2],sep=" --> ")))[,2],decreasing=T),]
write.table(epiicd9int,"Results/ICD9_comorbidities_categories.txt",quote=F,sep="\t",row.names=F,col.names=F)
## Save the different networks
complete<-network[which(as.numeric(network[,3])>=1),1:2]
write.table(complete,"Data/Network_analysis/PDN_complete.txt",quote=F,sep="\t",col.names=F,row.names=F)
pdnepi<-cbind(gsub("_.+","",pdn2$Epidemiological),gsub(".+_","",pdn2$Epidemiological))
write.table(pdnepi,"Data/Network_analysis/PDN_subset.txt",quote=F,sep="\t",col.names=F,row.names=F)
pdnpos<-cbind(gsub("_.+","",pdn2$Positive),gsub(".+_","",pdn2$Positive))
write.table(pdnpos,"Data/Network_analysis/ICD9_positive.txt",quote=F,sep="\t",col.names=F,row.names=F)
pdnneg<-cbind(gsub("_.+","",pdn2$Negative),gsub(".+_","",pdn2$Negative))
write.table(pdnneg,"Data/Network_analysis/ICD9_negative.txt",quote=F,sep="\t",col.names=F,row.names=F)

pdn<-pdn2$Epidemiological
save(pdn,file="Data/PDN_shared_diseases.Rdata")    

## @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ ##
## Comorbidity analyses for ICD9, ICD10, Soren and Barabasi ##
## @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ ##
load("Data/RR/ICD9_pRR.Rdata")
load("Data/RR/ICD9_nRR.Rdata")
positivos<-RRpos[which(as.numeric(RRpos[,5])>1),1:2]   ;   negativos<-RRneg[which(as.numeric(RRneg[,5])>1),1:2]
posipeg<-paste(positivos[,1],positivos[,2],sep="_")   ;   negapeg<-paste(negativos[,1],negativos[,2],sep="_")
posipeg2<-setdiff(posipeg,negapeg)   ;   negapeg2<-setdiff(negapeg,posipeg)   ;   posipeg<-posipeg2   ;   negapeg<-negapeg2
splited<-strsplit(posipeg,"_") ; pos<-c() ; for(a in 1:length(splited)){pos<-rbind(pos,c(splited[[a]][1],splited[[a]][2]))}
splited<-strsplit(negapeg,"_") ; neg<-c() ; for(a in 1:length(splited)){neg<-rbind(neg,c(splited[[a]][1],splited[[a]][2]))}
positive<-pos[-unique(c(grep("ICD9--",pos[,1]),grep("ICD9--",pos[,2]))),]
negative<-neg[-unique(c(grep("ICD9--",neg[,1]),grep("ICD9--",neg[,2]))),]
info<-read.csv2("Patients_information_table.txt",stringsAsFactors = F,sep="\t",header=T)
cambio<-info[,7]
names(cambio)<-gsub("_","-",info[,4])
positivo<-cbind(as.character(cambio[positive[,1]]),as.character(cambio[positive[,2]]))
negativo<-cbind(as.character(cambio[negative[,1]]),as.character(cambio[negative[,2]]))
ppos<-as.data.frame(table(paste(positivo[,1],positivo[,2],sep="_")))
splited<-strsplit(as.character(ppos[,1]),"_")
pinter<-c();for(a in 1:length(splited)){pinter<-rbind(pinter,c(splited[[a]][1],splited[[a]][2]))}
network<-cbind(pinter,as.character(ppos[,2]),"+")
pneg<-as.data.frame(table(paste(negativo[,1],negativo[,2],sep="_")))
splited<-strsplit(as.character(pneg[,1]),"_")
ninter<-c();for(a in 1:length(splited)){ninter<-rbind(ninter,c(splited[[a]][1],splited[[a]][2]))}
network<-rbind(network,cbind(ninter,as.character(pneg[,2]),"-"))
colnames(network)<-c("Interactor_1","Interactor_2","Intensity","Direction")
write.table(network,"Networks/ICD9_category_interactions.txt",quote=F,sep="\t",row.names=F)

## Write the table to represent the network of intra-disease patient-patient similarity network ##
## @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ ##
load("Data/Patient_patient_similarity_network.Rdata")
table(fisher_table[,3])
red<-fisher_table[which(fisher_table[,3]!=0),]
intra<-red[which(gsub("_.+","",as.character(red[,1]))==gsub("_.+","",as.character(red[,2]))),]
write.table(intra,"Networks/Patient_intra_disease_interactions.txt",quote=F,sep="\t",row.names=F)

## Number of overlaps and combinations detected in our analyses ##
## @@ @@ @@ @@ @@ @@ @@ @@ @@ @@  @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ ##
load("Data/Fisher_test_results.Rdata")
updown<-which(as.numeric(fisherresults$Fisher_values[,3])<=0.0001)
upup<-which(as.numeric(fisherresults$Fisher_values[,4])<=0.0001)
downup<-which(as.numeric(fisherresults$Fisher_values[,5])<=0.0001)
downdown<-which(as.numeric(fisherresults$Fisher_values[,6])<=0.0001)
install.packages("VennDiagram")
library("VennDiagram")
combinaciones<-list("UpUp"=upup,"DownUp"=downup,"UpDown"=updown,"DownDown"=downdown)
venn.diagram(combinaciones,"Plots/Number_of_overlaps_between_deg")

