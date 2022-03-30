# Patient-specific relative risk estimator
setwd("./")
load("Data/Patient_patient_similarity_network.Rdata") # fisher_table
red000005<-fisher_table  ;  enfermedad<-unique(gsub("_.+","",unique(c(red000005[,1],red000005[,2]))))
## Remove intra-disease interactions
uno<-gsub("_.+","",red000005[,1],fixed=F)  ;  dos<-gsub("_.+","",red000005[,2],fixed=F)  ;  intra<-which(uno==dos)  ;  inter000005<-red000005[-intra,]
tablon<-read.csv2("Patients_information_table.txt",stringsAsFactors = F,sep="\t",header=T)
pacientes<-tablon[,1]
## Patient identifiers
identidir<-c("Data/Identifiers/Patient_identifiers/")
identifiers<-list()
for(z in 1:length(list.files(identidir))){
  load(paste(identidir,list.files(identidir)[z],sep="")) ; for(q in 1:length(names(identificadores))){identifiers[[names(identificadores)[q]]]<-identificadores[[q]]}
}
load("Data/Identifiers/Disease.Rdata") # identificadores
for(a in 1:length(names(identifiers))){
  if(a!=3675){
    print(a)  ;  enf1<-identifiers[[a]]
    for(b in 1:length(enfermedad)){
      if(gsub("_.+","",names(identifiers)[a])!=enfermedad[b]){
        ## Declare objects
        nopresentesiendo<-c() ; presentesiendo<-c() ; presentenosiendo<-c() ; nopresentenosiendo<-c()
        nopresentesiendo1<-c() ; nopresentesiendo2<-c() ; nopresentenosiendo1<-c() ; nopresentenosiendo2<-c()
        ## Select the interactions of interest
        enf2<-identificadores[[enfermedad[b]]] ; par<-intersect(enf1,enf2) ; solo1<-setdiff(enf1,par) ; solo2<-setdiff(enf2,par)
        ## @ @ ##
        ## pRR ##
        siendo<-table(inter000005[par,3])
        # Not present being
        if(length(which(names(siendo)==0))>0){nopresentesiendo1<-as.numeric(siendo[which(names(siendo)==0)])}
        if(length(which(names(siendo)==(-1)))>0){nopresentesiendo2<-as.numeric(siendo[which(names(siendo)==(-1))])}
        if(length(nopresentesiendo1)>0){nopresentesiendo<-nopresentesiendo1}
        if(length(nopresentesiendo2)>0){nopresentesiendo<-nopresentesiendo2}
        if(length(nopresentesiendo1)>0 && length(nopresentesiendo2)>0){nopresentesiendo<-nopresentesiendo1+nopresentesiendo2}
        if(length(nopresentesiendo)==0){nopresentesiendo<-0}                                    ## NEW 22/11/2017
        # Present being
        if(length(which(names(siendo)==1))>0){presentesiendo<-as.numeric(siendo[which(names(siendo)==1)])}
        if(length(presentesiendo)==0){presentesiendo<-0}                                        ## NEW 22/11/2017
        ## Disease
        nosiendo<-table(inter000005[solo1,3])
        # No present not being
        if(length(which(names(nosiendo)==0))>0){nopresentenosiendo1<-as.numeric(nosiendo[which(names(nosiendo)==0)])}
        if(length(which(names(nosiendo)==(-1)))>0){nopresentenosiendo2<-as.numeric(nosiendo[which(names(nosiendo)==(-1))])}
        if(length(nopresentenosiendo1)>0){nopresentenosiendo<-nopresentenosiendo1}
        if(length(nopresentenosiendo2)>0){nopresentenosiendo<-nopresentenosiendo2}
        if(length(nopresentenosiendo1)>0 && length(nopresentenosiendo2)>0){nopresentenosiendo<-nopresentenosiendo1+nopresentenosiendo2}
        if(length(nopresentenosiendo)==0){nopresentenosiendo<-0}                                    ## NEW 22/11/2017
        # Present not being
        if(length(which(names(nosiendo)==1))>0){presentenosiendo<-as.numeric(nosiendo[which(names(nosiendo)==1)])}
        if(length(presentenosiendo)==0){presentenosiendo<-0}                                    ## NEW 22/11/2017
        peventwhenexposed<-presentesiendo/(presentesiendo+nopresentesiendo)
        peventwhennotexposed<-presentenosiendo/(presentenosiendo+nopresentenosiendo)
        division<-peventwhenexposed/peventwhennotexposed
        if(length(division)>0){
          raiz<-sqrt(((nopresentesiendo/presentesiendo)/(nopresentesiendo+presentesiendo))+((nopresentenosiendo/presentenosiendo)/(nopresentenosiendo+presentenosiendo)))
          sup<-exp(log(division)+(1.96*raiz))  ;  inf<-exp(log(division)-(1.96*raiz))
          # RRpos<-rbind(RRpos,c(names(identifiers)[a],enfermedad[b],division,sup,inf))
          jun<-c(names(identifiers)[a],enfermedad[b],division,sup,inf)
          jun<-t(jun)
          write.table(jun,"Data/RR/Patient_pRR.txt",quote=F,sep="\t",row.names=F,col.names=F,append = T)
        }
        ## @ @ ##
        ## nRR ##
        ## Declare objectes
        nopresentesiendo<-c() ; presentesiendo<-c() ; presentenosiendo<-c() ; nopresentenosiendo<-c()
        nopresentesiendo1<-c() ; nopresentesiendo2<-c() ; nopresentenosiendo1<-c() ; nopresentenosiendo2<-c()
        siendo<-table(inter000005[par,3])
        # Not present being
        if(length(which(names(siendo)==0))>0){nopresentesiendo1<-as.numeric(siendo[which(names(siendo)==0)])}
        if(length(which(names(siendo)==1))>0){nopresentesiendo2<-as.numeric(siendo[which(names(siendo)==1)])}
        if(length(nopresentesiendo1)>0){nopresentesiendo<-nopresentesiendo1}
        if(length(nopresentesiendo2)>0){nopresentesiendo<-nopresentesiendo2}
        if(length(nopresentesiendo1)>0 && length(nopresentesiendo2)>0){nopresentesiendo<-nopresentesiendo1+nopresentesiendo2}
        if(length(nopresentesiendo)==0){nopresentesiendo<-0}                                    ## NEW 22/11/2017
        # Present being
        if(length(which(names(siendo)==(-1)))>0){presentesiendo<-as.numeric(siendo[which(names(siendo)==(-1))])}
        if(length(presentesiendo)==0){presentesiendo<-0}                                    ## NEW 22/11/2017
        ## Disease
        nosiendo<-table(inter000005[solo1,3])
        # Not present not being
        if(length(which(names(nosiendo)==0))>0){nopresentenosiendo1<-as.numeric(nosiendo[which(names(nosiendo)==0)])}
        if(length(which(names(nosiendo)==1))>0){nopresentenosiendo2<-as.numeric(nosiendo[which(names(nosiendo)==1)])}
        if(length(nopresentenosiendo1)>0){nopresentenosiendo<-nopresentenosiendo1}
        if(length(nopresentenosiendo2)>0){nopresentenosiendo<-nopresentenosiendo2}
        if(length(nopresentenosiendo1)>0 && length(nopresentenosiendo2)>0){nopresentenosiendo<-nopresentenosiendo1+nopresentenosiendo2}
        if(length(nopresentenosiendo)==0){nopresentenosiendo<-0}                                    ## NEW 22/11/2017
        # Present not being
        if(length(which(names(nosiendo)==(-1)))>0){presentenosiendo<-as.numeric(nosiendo[which(names(nosiendo)==(-1))])}
        if(length(presentenosiendo)==0){presentenosiendo<-0}                                    ## NEW 22/11/2017
        peventwhenexposed<-presentesiendo/(presentesiendo+nopresentesiendo)
        peventwhennotexposed<-presentenosiendo/(presentenosiendo+nopresentenosiendo)
        division<-peventwhenexposed/peventwhennotexposed
        if(length(division)>0){
          raiz<-sqrt(((nopresentesiendo/presentesiendo)/(nopresentesiendo+presentesiendo))+((nopresentenosiendo/presentenosiendo)/(nopresentenosiendo+presentenosiendo)))
          sup<-exp(log(division)+(1.96*raiz))  ;  inf<-exp(log(division)-(1.96*raiz))
          # RRneg<-rbind(RRneg,c(names(identifiers)[a],enfermedad[b],division,sup,inf))
          jun<-c(names(identifiers)[a],enfermedad[b],division,sup,inf)
          jun<-t(jun)
          write.table(jun,"Data/RR/Patient_nRR.txt",quote=F,sep="\t",row.names=F,col.names=F,append = T)
        }
      }
    }
  }
}
print("FIN!!")

RRpos<-read.csv2("Data/RR/Patient_pRR.txt",stringsAsFactors = F,sep="\t",header=F) ; RRneg<-read.csv2("Data/RR/Patient_nRR.txt",stringsAsFactors = F,sep="\t",header=F)
RRspositivas<-RRpos[which(as.numeric(RRpos[,5])>1),] ; RRsnegativas<-RRneg[which(as.numeric(RRneg[,5])>1),]
colnames(RRsnegativas)<-c("Patient","Disease","RR","CI_95%_up","CI_95%_down") ; colnames(RRspositivas)<-c("Patient","Disease","RR","CI_95%_up","CI_95%_down")

neg<-RRsnegativas[,1:3] ; pos<-RRspositivas[,1:3]
pasneg<-paste(neg[,1],neg[,2],sep="_") ; paspos<-paste(pos[,1],pos[,2],sep="_")
quit<-intersect(pasneg,paspos)
quitpos<-c();for(a in 1:length(quit)){quitpos<-c(quitpos,which(paspos==quit[a]))}  ;  quitneg<-c();for(a in 1:length(quit)){quitneg<-c(quitneg,which(pasneg==quit[a]))}
negativos<-neg[-quitneg,] ; positivos<-pos[-quitpos,]
po<-cbind(positivos,1) ; colnames(po)[4]<-"Interaction" ; ne<-cbind(negativos,-1) ; colnames(ne)[4]<-"Interaction"
rr1<-rbind(po,ne)

## Alzheimer's disease patients protected against NSCLC ##
## @@ @@ @@ @@ @@ @@ @@ @@ @ @@ @@ @@ @@ @@ @@ @@ @@ @@ ##
adlc<-rr1[intersect(grep("Alzheimer",rr1[,1]),grep("NSCLC",rr1[,2])),]
tabla<-read.csv2("Patients_information_table.txt",stringsAsFactors = F,sep="\t",header=T)
enfermedades<-unique(tabla$Disease)
sub1<-tabla$Cluster_info_cluster_specific ; names(sub1)<-tabla$Patient
dispat<-"NSCLC"
disint<-"AlzheimersDisease"
load("Data/Top_500_sDEGs_names.Rdata")
load("Data/Drugs.Rdata")
direct<-list() ; inverse<-list()
for(a in 1:length(names(drugs))){
  direct[[names(drugs)[a]]]<-drugs[[a]]$Similar
  inverse[[names(drugs)[a]]]<-drugs[[a]]$Opposite
}
relevantes<-function(dispat,disint,rr=rr1,sub=sub1,de=degs,direct1=direct,inverse1=inverse){
  if(length(intersect(grep(paste("^",dispat,"_",sep=""),rr[,1],ignore.case = T),grep(paste("^",disint,"$",sep=""),rr[,2],ignore.case = T)))>0){
    irr<-rr[intersect(grep(paste("^",dispat,"_",sep=""),rr[,1],ignore.case = T),grep(paste("^",disint,"$",sep=""),rr[,2],ignore.case = T)),]
    irrp<-c();irrn<-c()
    if(length(irr)>4){
      if(length(which(irr[,4]==1))>0){irrp<-irr[which(irr[,4]==1),] ; if(length(irrp)>4){irrp<-irrp[order(as.numeric(irrp[,3]),decreasing = T),]}}
      if(length(which(irr[,4]==-1))>0){irrn<-irr[which(irr[,4]==-1),] ; if(length(irrn)>4){irrn<-irrn[order(as.numeric(irrn[,3]),decreasing = T),]}}}
    if(length(irr)==4){
      if(length(which(irr[4]==1))>0){irrp<-irr ; if(length(irrp)>4){irrp<-irrp[order(as.numeric(irrp[,3]),decreasing = T),]}}
      if(length(which(irr[4]==-1))>0){irrn<-irr ; if(length(irrn)>4){irrn<-irrn[order(as.numeric(irrn[,3]),decreasing = T),]}}
    }
    irr<-rbind(irrp,irrn)
    irr<-cbind(as.character(sub[irr[,1]]),irr)
    colnames(irr)[1]<-"Patients_Subgroup"
    posrisk<-paste(round((length(which(irr[,5]==1))/length(grep(paste("^",dispat,"_",sep=""),names(sub))))*100,3),"%",sep="")
    negrisk<-paste(round((length(which(irr[,5]==-1))/length(grep(paste("^",dispat,"_",sep=""),names(sub))))*100,3),"%",sep="")
    ## Shared genes
    top10_prr<-list()  ;  drug_top10_prr<-list()
    if(round((length(which(irr[,5]==1))/length(grep(dispat,sub)))*100,3)>0){
      selpat<-irr[which(irr[,5]==1),2]
      ## genes
      upl<-list();downl<-list() ; ups<-c();downs<-c()
      for(a in 1:length(selpat)){upl[[a]]<-de[[selpat[a]]]$up  ;  ups<-c(ups,de[[selpat[a]]]$up)  ;  downl[[a]]<-de[[selpat[a]]]$down  ;  downs<-c(downs,de[[selpat[a]]]$down)}
      comup<-Reduce(intersect,upl)
      comdo<-Reduce(intersect,downl)
      up<-table(ups)[order(as.numeric(table(ups)),decreasing=T)][1:10]/length(selpat)
      down<-table(downs)[order(as.numeric(table(downs)),decreasing=T)][1:10]/length(selpat)
      top10_prr$up<-up ; top10_prr$down<-down
      ## drugs
      upl<-list();downl<-list()  ;  ups<-c();downs<-c()
      for(a in 1:length(selpat)){upl[[a]]<-direct1[[selpat[a]]]  ;  ups<-c(ups,direct1[[selpat[a]]])  ;  downl[[a]]<-inverse1[[selpat[a]]]  ;  downs<-c(downs,inverse1[[selpat[a]]])}
      comup<-Reduce(intersect,upl)
      comdo<-Reduce(intersect,downl)
      up<-table(ups)[order(as.numeric(table(ups)),decreasing=T)][1:10]/length(selpat)
      down<-table(downs)[order(as.numeric(table(downs)),decreasing=T)][1:10]/length(selpat)
      drug_top10_prr$direct<-up ; drug_top10_prr$inverse<-down
    }
    top10_nrr<-list()  ;  drug_top10_nrr<-list()
    if(round((length(which(irr[,5]==-1))/length(grep(dispat,sub)))*100,3)>0){
      selpat<-irr[which(irr[,5]==-1),2]
      upl<-list();downl<-list()
      ups<-c();downs<-c()
      for(a in 1:length(selpat)){upl[[a]]<-de[[selpat[a]]]$up  ;  ups<-c(ups,de[[selpat[a]]]$up)  ;  downl[[a]]<-de[[selpat[a]]]$down  ;  downs<-c(downs,de[[selpat[a]]]$down)}
      comup<-Reduce(intersect,upl)
      comdo<-Reduce(intersect,downl)
      up<-table(ups)[order(as.numeric(table(ups)),decreasing=T)][1:10]/length(selpat)
      down<-table(downs)[order(as.numeric(table(downs)),decreasing=T)][1:10]/length(selpat)
      top10_nrr$up<-up
      top10_nrr$down<-down
      ## drugs
      upl<-list();downl<-list()  ;  ups<-c();downs<-c()
      for(a in 1:length(selpat)){upl[[a]]<-direct1[[selpat[a]]]  ;  ups<-c(ups,direct1[[selpat[a]]])  ;  downl[[a]]<-inverse1[[selpat[a]]]  ;  downs<-c(downs,inverse1[[selpat[a]]])}
      comup<-Reduce(intersect,upl)
      comdo<-Reduce(intersect,downl)
      up<-table(ups)[order(as.numeric(table(ups)),decreasing=T)][1:10]/length(selpat)
      down<-table(downs)[order(as.numeric(table(downs)),decreasing=T)][1:10]/length(selpat)
      drug_top10_nrr$direct<-up
      drug_top10_nrr$inverse<-down
    }
    if(round((length(which(irr[,5]==1))/length(grep(dispat,sub)))*100,3)>0 && round((length(which(irr[,5]==-1))/length(grep(dispat,sub)))*100,3)>0){
      resultados<-list("Interactions"=irr,"Positive_risk"=posrisk,"Negative_risk"=negrisk,"top10_nrr"=top10_nrr,"top10_prr"=top10_prr,
                       "drug_top10_prr"=drug_top10_prr,"drug_top10_nrr"=drug_top10_nrr)
    }
    if(round((length(which(irr[,5]==1))/length(grep(dispat,sub)))*100,3)>0 && round((length(which(irr[,5]==-1))/length(grep(dispat,sub)))*100,3)==0){
      resultados<-list("Interactions"=irr,"Positive_risk"=posrisk,"Negative_risk"=negrisk,"top10_prr"=top10_prr,"drug_top10_prr"=drug_top10_prr)
    }
    if(round((length(which(irr[,5]==1))/length(grep(dispat,sub)))*100,3)==0 && round((length(which(irr[,5]==-1))/length(grep(dispat,sub)))*100,3)>0){
      resultados<-list("Interactions"=irr,"Positive_risk"=posrisk,"Negative_risk"=negrisk,"top10_nrr"=top10_nrr,"drug_top10_nrr"=drug_top10_nrr)
    }
    return(resultados)
  }
}

an<-relevantes("AlzheimersDisease","NSCLC")
na<-relevantes("NSCLC","AlzheimersDisease")


## Extract the drug pool ##
completedrugs<-c()  ;   for(a in 1:length(names(direct))){completedrugs<-c(completedrugs,direct[[a]],inverse[[a]])}   ;   completedrugs<-unique(completedrugs)
## Obtain a list of ranked drugs and directions for each disease with the percentage of patients taking them
disdrug<-list()
for(a in 1:length(enfermedades)){
  # a<-2
  paterest<-tabla[which(tabla$Disease==enfermedades[a]),1]
  drojas<-c()
  for(b in 1:length(paterest)){
    # b<-1
    drojas<-c(drojas,paste(direct[[paterest[b]]],"__1",sep=""),paste(direct[[paterest[b]]],"__-1",sep=""))
  }
  sun<-table(drojas)[order(as.numeric(table(drojas)),decreasing=T)]/length(paterest)
  sun2<-as.data.frame(sun)
  sunp<-sun2[grep("__1",as.character(sun2[,1])),]
  sunn<-sun2[grep("__-1",as.character(sun2[,1])),]
  comors<-rbind(cbind(gsub("__1","",sunp[,1]),sunp[,2],1),cbind(gsub("__1","",sunn[,1]),sunn[,2],-1))
  comors<-comors[order(as.numeric(comors[,2]),decreasing=T),]
  colnames(comors)<-c("Drugs","Presence","Direction")
  disdrug[[enfermedades[a]]]<-comors
}

## Drug finder (if it exists in LINCS)
completedrugs[grep("galantamine",completedrugs,ignore.case = T)]

## Alzheimer drugs 
addrugs<-c("galantamine","risperidone","rivastigmine","donepezil","ziprasidone")

## Look for patient-specific comorbidities
alzheimerspat<-tabla[grep("Alzheimer",tabla[,1]),1]
buscapat<-function(patient,value=direct,threshold=0.6,rr=rr1){
  if(length(which(rr[,1]==patient))>0){
    son<-rr[which(rr[,1]==patient),]
    sonp<-c()  ; sonn<-c()
    if(sum(dim(son))>4){
      if(length(which(son[,4]==1))>0){
        if(length(which(son[,4]==1))>1){sonp<-son[which(son[,4]==1),] ; sonp<-sonp[order(as.numeric(sonp[,4]),decreasing = T),]}
        if(length(which(son[,4]==1))==1){sonp<-son[which(son[,4]==1),]}
      }
      if(length(which(son[,4]==-1))>0){
        if(length(which(son[,4]==-1))>1){sonn<-son[which(son[,4]==-1),] ; sonn<-sonn[order(as.numeric(sonn[,4]),decreasing = T),]}
        if(length(which(son[,4]==-1))==1){sonn<-son[which(son[,4]==-1),]}
      }
      tabpat<-as.data.frame(rbind(sonp,sonn))
      enfgas<-list()
      # Buscamos los casos de comorbilidad directa
      if(length(which(tabpat[,4]==1))>0){
        # Las seleccionamos
        tapi<-tabpat[which(tabpat[,4]==1),]
        for(a in 1:length(tapi[,1])){
          # Drogas de la enfermedad 
          todas<-disdrug[which(names(disdrug)==as.character(tapi[a,2]))][[1]]
          # Cogemos las que causan la segunda enfermedad
          todas<-todas[which(todas[,3]==1),]
          # Al menos que esten en un x% de los pacientes en dicha direccion
          cogidas<-todas[which(as.numeric(todas[,2])>=threshold),]
          tenemos<-c()
          # Solapa alguna con la de nuestros pacientes?
          # Dependiendo de value miraremos las directas o inversas del paciente
          if(sum(dim(cogidas))>3){rownames(cogidas)<-cogidas[,1]  ;  tenemos<-intersect(value[[patient]],cogidas[,1])}
          if(sum(dim(cogidas))==3){tenemos<-intersect(value[[patient]],cogidas[1])}
          # Tenemos mas de un farmaco en comun
          if(length(tenemos)>1){puts<-cogidas[tenemos,][order(as.numeric(cogidas[tenemos,2]),decreasing = T),] ; enfgas[[as.character(tapi[a,2])]]<-puts}
          # Tenemos un solo farmaco en comun
          if(length(tenemos)==1){
            if(length(cogidas)>3){puts<-cogidas[tenemos,] ; enfgas[[as.character(tapi[a,2])]]<-puts}
            if(length(cogidas)==3){puts<-cogidas ; enfgas[[as.character(tapi[a,2])]]<-puts}
          }
        }
      }
      resultados<-list("Tabla"=tabpat,"Drugs"=enfgas)
      return(resultados)
    }
  }
  if(length(which(rr[,1]==patient))==0){print("No tiene comorbilidades")}
}

## Look for patients with drugs of interest to be avoided
toditas<-c()
for(a in 1:length(alzheimerspat)){
  # a<-2
  res<-buscapat(alzheimerspat[a],inverse) ## change this if we want to see a different association
  if(class(res)!="character"){
    if(length(res$Drugs)>0){
      drogas<-res$Drugs
      ndr<-names(drogas)
      td<-c()
      for(b in 1:length(ndr)){
        # b<-1
        if(length(drogas[[ndr[b]]])==3){td<-rbind(td,c(ndr[b],drogas[[ndr[b]]]))}
        if(length(drogas[[ndr[b]]])>3){td<-rbind(td,cbind(ndr[b],drogas[[ndr[b]]]))}
      }
      toditas<-rbind(toditas,cbind(alzheimerspat[a],td))
    }
  }
}
completedrugs[grep("olanzapine",completedrugs,ignore.case = T)]

toditas[which(as.character(toditas[,3])=="cyproterone"),]
buscapat("AlzheimersDisease_GSE5281_Patient_8",inverse)


