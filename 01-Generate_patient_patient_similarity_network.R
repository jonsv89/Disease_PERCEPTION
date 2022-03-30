## CODE FOR THE DIFFERENTIAL EXPRESSION ANALYSIS OF EACH CASE SAMPLE AGAINST ALL THE CONTROL SAMPLES FROM THE SAME STUDY ##
## Developed by Jon Sanchez-Valle
## Barcelona Supercomputing Center
## Life Science Department
## Computational Biology Group
## Email: jon.sanchez@bsc.es

## load the essential packages 
library("affy");library(frma);library("limma");library("genefilter");require(graphics);library("hgu133plus2frmavecs");library("massiR");library("hgu133plus2.db")
library("cluster");library("NbClust");library(e1071)
## generate required directories
if ("ToptableIQR"%in%list.files("Data/") == FALSE){dir.create("Data/ToptableIQR")}
if ("ExpressionSet_Together"%in%list.files("Data/") == FALSE){dir.create("Data/ExpressionSet_Together")}
if ("ExpressionGenes_Together"%in%list.files("Data/") == FALSE){dir.create("Data/ExpressionGenes_Together")}
if ("Correspondences"%in%list.files("Data/") == FALSE){dir.create("Data/Correspondences")}
if ("PatientIQR"%in%list.files("Data/") == FALSE){dir.create("Data/PatientIQR")}
if ("Patienttvalues"%in%list.files("Data/") == FALSE){dir.create("Data/Patienttvalues")}
if ("Networks"%in%list.files() == FALSE){dir.create("Networks")}
if ("Data/ExpressionSetfiles"%in%list.files("Data/") == FALSE){dir.create("Data/ExpressionSetfiles")}


## set seed ##
## @@ @@ @@ ##
set.seed(1)

                       ## @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@  @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ ##
                       ##                       Generate gene expression matrices for each disease separated by studies                      ##
                       ## conduct differential expression analyses comparing in each study all the case samples with all the control samples ##
                       ## @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@  @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ ##

rawdatadir<-"./Data/Raw_data/"
rawdatafiles<-list.files(rawdatadir)
## load gene symbol - probe association and the probes that must be removed 
load("Data/Remove_Symbols_u133plus2.Rdata") # remo
load("Data/Gene_Symbols_u133plus2.Rdata") # simbolos
## indicate which diseases have been already analyzed 
if("Analyzed_diseases.txt"%in%list.files("Data/") == TRUE){echos<-read.table("Data/Analyzed_diseases.txt",stringsAsFactors = F)[,1];rawdatafiles<-setdiff(rawdatafiles,echos);dones<-echos}
if("Analyzed_diseases.txt"%in%list.files("Data/") == FALSE){dones<-c()}
print("We start the loop!")
for (a in 1:length(rawdatafiles)){
  inicio<-Sys.time()
  ## print information regarding the dataset that is going to be analyzed
  print(paste(rawdatafiles[a],"Control:",length(list.files(paste("./",rawdatadir,"/",rawdatafiles[a],"/Control/",sep=""))),"Case:",
              length(list.files(paste("./",rawdatadir,"/",rawdatafiles[a],"/Case/",sep=""))),sep="      "))
  control<-list.files(paste("./",rawdatadir,"/",rawdatafiles[a],"/Control",sep=""))
  path_con<-paste(rawdatadir,rawdatafiles[a],"/Control/",control,sep="")
  case<-list.files(paste("./",rawdatadir,"/",rawdatafiles[a],"/Case",sep=""))
  path_ca<-paste(rawdatadir,rawdatafiles[a],"/Case/",case,sep="")
  path_read<-c(path_con, path_ca)
  controls<-cbind(control,rep("control",length(control)))
  cases<-cbind(case,rep("case",length(case)))
  targets<-rbind(controls, cases)
  affyBatch<-ReadAffy(filenames=path_read)
  expSet<-frma(affyBatch)
  save(expSet,file=paste("Data/ExpressionSetfiles/",rawdatafiles[a],".Rdata",sep=""))
  barcoded<-barcode(expSet)
  save(barcoded,file=paste("Data/Barcodes/",rawdatafiles[a],".Rdata",sep=""))
  comoda<-t(as.data.frame(expSet))
  ## save the normalized expression matrix (with probes)
  save(comoda,file=paste("Data/ExpressionSet_Together/",rawdatafiles[a],".Rdata",sep=""))
  comoda1<-comoda[-length(comoda[,1]),]
  ## remove probes not associated to gene symbols
  comoda2<-comoda1[-remo,]
  coluna<-colnames(comoda2)
  genes<-unique(simbolos)
  expression<-c()
  ## convert the probes' expression matrix into a gene expression matrix, median values are calculated when seveeral probes refer to the same gene symbol
  for (c in 1:length(genes)){
    com<-which(simbolos==as.character(genes[c]))
    if (length(com)>1){
      iq<-apply(comoda[com,],2,median)
    }
    if (length(com)==1){
      iq<-comoda[com,]
    }
    expression<-rbind(expression,iq)
  }
  colnames(expression)<-coluna
  rownames(expression)<-genes
  columna<-c(control,case)
  ## change samples names
  newpatientnames<-c(paste("Control_",1:length(control),sep=""),paste("Patient_",1:length(case),sep=""))
  correspondencia<-cbind(columna,newpatientnames)
  colnames(correspondencia)<-c("Original_name","Patient_name")
  ## save patient name - sample ID information
  save(correspondencia,file=paste("Data/Correspondences/",rawdatafiles[a],".Rdata",sep=""))
  colnames(expression)<-newpatientnames
  controls<-cbind(paste("Control_",1:length(control),sep=""),rep("control",length(control)))
  cases<-cbind(paste("Patient_",1:length(case),sep=""),rep("case",length(case)))
  targets<-rbind(controls, cases)
  ## save gene expression matrix with the new names
  save(expression,file=paste("Data/ExpressionGenes_Together/",rawdatafiles[a],".Rdata",sep=""))
  ## conduct differential gene expression analysis using LIMMA
  design<-cbind(CONTROL=c(rep(1,length(control)),rep(0,length(case))), CASE=c(rep(0,length(control)),rep(1,length(case))))
  rownames(design)<-targets[,1]
  cont.matrix<-makeContrasts(CASEvsCONTROL=CASE-CONTROL,levels=design)
  ## get DEGs from IQR
  fit<-lmFit(expression,design)
  fit2<-contrasts.fit(fit, cont.matrix)
  fit2<-eBayes(fit2)
  toptableIQR<-topTable(fit2, number=length(fit$Amean), adjust.method="BH", sort.by="p")
  ## save all patient vs. all control differential expression results
  write.table(toptableIQR,paste("./Data/ToptableIQR/",rawdatafiles[a],".txt",sep=""),sep="\t",quote=F)
  print(paste(a,"  out of  ",length(rawdatafiles),":  ",rawdatafiles[a],sep=""))
  dones<-c(dones,rawdatafiles[a])
  ## save the name of the analyzed disease (with its' corresponding study)
  write.table(dones,"Data/Analyzed_diseases.txt",quote=F,col.names=F,row.names=F,sep="\t")
  final<-Sys.time()
  diferencia<-final-inicio
  print(diferencia)
}
print("Gene-expression matrices and global differential expression analyses finished!!")


                        ## @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ ##
                        ## Conduct differential expression analyses comparing each case sample with all the control samples from the same study ##
                        ## @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ ##

expresiondir<-c("Data/ExpressionGenes_Together/")
expresionfiles<-list.files(expresiondir)
nombrecitos<-gsub(".Rdata","",expresionfiles,fixed = T)
## avoid re-analyzing patients
if(length(list.files("Data/Patienttvalues/"))>0){
  pacientes<-list.files("Data/Patienttvalues/")
  if(length(pacientes)>0){
    hecho<-c()
    dones<-gsub("_Patient.+",".Rdata",pacientes)
    done<-unique(dones)
    for(a in done){
      load(paste(expresiondir,a,sep=""))
      if(length(which(dones==a))==length(grep("Patient",colnames(expression)))){
        hecho<-c(hecho,a)
      }
    }
  }
  if(length(pacientes)==0){
    hecho<-c()
  }
  expresionfiles<-setdiff(expresionfiles,hecho)
}
## start patient-specific differential expression analysis
for(a in 1:length(expresionfiles)){
  load(paste(expresiondir,expresionfiles[a],sep=""))
  patientes<-grep("Patient",colnames(expression))
  controles<-grep("Control",colnames(expression))
  for(b in 1:length(patientes)){
    expresion<-expression[,c(controles,patientes[b])]
    control<-colnames(expression)[controles]
    case<-colnames(expression)[patientes[b]]
    controls<-cbind(paste("Control_",1:length(control),sep=""),rep("control",length(control)))
    cases<-cbind(paste("Patient_",1:length(case),sep=""),rep("case",length(case)))
    targets<-rbind(controls, cases)
    design<-cbind(CONTROL=c(rep(1,length(control)),rep(0,length(case))), CASE=c(rep(0,length(control)),rep(1,length(case))))
    rownames(design)<-targets[,1]
    cont.matrix<-makeContrasts(CASEvsCONTROL=CASE-CONTROL,levels=design)
    fit<-lmFit(expresion,design)  ##getting DEGs from IQR
    fit2<-contrasts.fit(fit, cont.matrix)
    fit2<-eBayes(fit2)
    toptableIQR<-topTable(fit2, number=length(fit$Amean), adjust.method="BH", sort.by="p")
    write.table(toptableIQR,paste("./Data/PatientIQR/",nombrecitos[a],"_",case,".txt",sep=""),sep="\t",quote=F)
    tvalues<-toptableIQR[,3]
    names(tvalues)<-rownames(toptableIQR)
    save(tvalues,file=paste("./Data/Patienttvalues/",nombrecitos[a],"_",case,".Rdata",sep=""))
  }
  print(paste(expresionfiles[a],"  we have finished  ",a,"  out of  ",length(expresionfiles),sep=""))
}
print("Patient-specific differential gene expression finished!!")

## Extract top X up- and down-regulated genes and generate a matrix of tvalues for LINCS analysis ##
## @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ ##
## duration: 35.88323 mins
topini<-Sys.time()
directory<-"Data/Patienttvalues/"
## indicate the number of DEG we want to select
seleccion<-500
patients<-list.files(directory)
nombres<-gsub(".Rdata","",patients,fixed=T)
## avoid re-analyzing patients
if (paste("Top_",seleccion,"_sDEGs.Rdata",sep="")%in%list.files("Data/") == TRUE){
  load(paste("Data/Top_",seleccion,"_sDEGs.Rdata",sep=""))
  hechos<-paste(names(sdegs),".Rdata",sep="")
  echos<-names(sdegs)
  patients<-setdiff(patients,hechos)
  nombres<-setdiff(nombres,echos) ; print("We do have sDEGs")
}
if (paste("Top_",seleccion,"_sDEGs.Rdata",sep="")%in%list.files("Data/") == FALSE){sdegs<-list() ; print("We do not have sDEGs")}
## select the up- and down-regulated genes
if(length(patients)>0){
  load(paste(directory,patients[1],sep=""))
  genes<-names(tvalues)
  tablatvalores<-c()
  for(a in 1:length(patients)){
    # a<-1
    load(paste(directory,patients[a],sep=""))
    tablatvalores<-cbind(tablatvalores,as.numeric(tvalues[genes]))
    tvalores<-sort(tvalues)
    sdegs[[nombres[a]]]$down<-tvalores[1:seleccion]
    sdegs[[nombres[a]]]$up<-tvalores[(length(tvalores)-(seleccion-1)):length(tvalores)]
    print(paste(a,"  out of  ",length(patients),"  patients",sep=""))
  }
  rownames(tablatvalores)<-genes
  colnames(tablatvalores)<-nombres
  save(sdegs,file=paste("Data/Top_",seleccion,"_sDEGs.Rdata",sep=""))
  saveRDS(tablatvalores,file="Data/LINCS_tvalues_patients.rds")
  print(paste("Top ",seleccion,"genes selected!!",sep="  "))
}
topfin<-Sys.time()
print(topfin-topini)

## Generate a list with DEGs names ##
## @@ @@ @@ @@ @@ @ @@ @@ @@ @@ @@ ##
## duration: 1.329571 mins
indeg<-Sys.time()
patients<-list.files(directory)
nombres<-gsub(".Rdata","",patients,fixed=T)
## avoid re-analyzing patients
if (paste("Top_",seleccion,"_sDEGs_names.Rdata",sep="")%in%list.files("Data/") == FALSE){degs<-list() ; print("We do not have the list with sDEGs names")}
if (paste("Top_",seleccion,"_sDEGs_names.Rdata",sep="")%in%list.files("Data/") == TRUE){
  load(paste("Data/Top_",seleccion,"_sDEGs_names.Rdata",sep=""))
  hechos<-paste(names(degs),".Rdata",sep="")
  echos<-names(degs)
  patients<-setdiff(patients,hechos)
  nombres<-setdiff(nombres,echos)
  print("We have the list with sDEGs names")
}
## select the names of the up- and down-regulated genes
if(length(patients)>0){
  for(a in 1:length(patients)){
    load(paste(directory,patients[a],sep=""))
    tvalores<-sort(tvalues)
    degs[[nombres[a]]]$down<-names(tvalores[1:seleccion])
    degs[[nombres[a]]]$up<-names(tvalores[(length(tvalores)-(seleccion-1)):length(tvalores)])
    print(paste(a,"  out of  ",length(patients),"  finished",sep=""))
  }
  save(degs,file=paste("Data/Top_",seleccion,"_sDEGs_names.Rdata",sep=""))
}
fideg<-Sys.time()
print(fideg-indeg)
fisin<-Sys.time()
load(paste("Data/Top_",seleccion,"_sDEGs_names.Rdata",sep=""))

## Fisher test analysis ##
## @@ @@ @@ @@ @@ @@ @@ ##
## duration: 17.61 hours
background<-20502
load("Data/Top_500_sDEGs_names.Rdata")
fisin<-Sys.time()
todos<-c()
toditos<-c()
for (a in 1:(length(names(degs))-1)){
  positi1<-degs[[a]]$up ; negati1<-degs[[a]]$down ; pega<-c() ; pegados<-c()
  for (b in (a+1):length(names(degs))){
    positi2<-degs[[b]]$up ; negati2<-degs[[b]]$down
    interspd <- length(intersect(positi1,positi2))
    kkpd <- matrix(c(interspd,length(positi1)-interspd,length(positi2)-interspd,background+interspd-length(positi1)-length(positi2)),nrow=2,ncol=2)
    fispd<-fisher.test(kkpd,alternative="greater") $p.value
    intersnd <- length(intersect(negati1,negati2))
    kknd <- matrix(c(intersnd,length(negati1)-intersnd,length(negati2)-intersnd,background+intersnd-length(negati1)-length(negati2)),nrow=2,ncol=2)
    fisnd<-fisher.test(kknd,alternative="greater") $p.value
    interspi <- length(intersect(positi1,negati2))
    kkpi <- matrix(c(interspi,length(positi1)-interspi,length(negati2)-interspi,background+interspi-length(positi1)-length(negati2)),nrow=2,ncol=2)
    fispi<-fisher.test(kkpi,alternative="greater") $p.value
    intersni <- length(intersect(negati1,positi2))
    kkni <- matrix(c(intersni,length(negati1)-intersni,length(positi2)-intersni,background+intersni-length(negati1)-length(positi2)),nrow=2,ncol=2)
    fisni<-fisher.test(kkni,alternative="greater") $p.value
    junt<-c(names(degs)[a],names(degs)[b],fispi,fispd,fisni,fisnd)
    jun<-t(junt)
    write.table(jun,"Fisher_overlaps.txt",quote=F,sep="\t",row.names=F,col.names=F,append = T)
  }
  print(a/(length(names(degs))-1))
}
print("Finished!!") ; fisfin<-Sys.time() ; print(fisfin-fisin)
## read the file and generate the network
fis<-read.csv("Fisher_overlaps.txt",stringsAsFactors = F,sep="\t",header=F)
threshold<-0.001
negs<-list("Pos_Neg"=which(fis[,3]<=threshold),"Neg_Pos"=which(fis[,5]<=threshold),"Pos_Pos"=which(fis[,4]>threshold),"Neg_Neg"=which(fis[,6]>threshold))
poss<-list("Pos_Neg"=which(fis[,3]>threshold),"Neg_Pos"=which(fis[,5]>threshold),"Pos_Pos"=which(fis[,4]<=threshold),"Neg_Neg"=which(fis[,6]<=threshold))
interactions<-rep(0,length(fis[,1])) ; interactions[Reduce(intersect,negs)]<--1 ; interactions[Reduce(intersect,poss)]<-1
## save the patient-patient similarity network including 0s
fisher_table<-cbind(fis[,1:2],interactions) ; colnames(fisher_table)<-c("Disease_1","Disease_2","Interaction")
save(fisher_table,file="Data/Patient_patient_similarity_network.Rdata")
network<-fisher_table[which(fisher_table[,3]!=0),]
## save the generated "patient-patient similarity network"
write.table(network,"Networks/Patient_patient_similarity_network.txt",sep="\t",quote=F,row.names=F)

## Generate a matrix of 1s, 0s and -1s with patients in rows and genes in columns ##
## @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @ @ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ ##
## estimated duration: 2.74 hours
starting<-Sys.time()
patients<-list.files("Data/Patienttvalues/") ; nome<-gsub(".Rdata","",patients,fixed=T) ; onesandzerosm<-list()
for(a in 1:length(patients)){
  load(paste("Data/Patienttvalues/",patients[a],sep=""))
  updown<-rep(0,length(tvalues)) ; updown[1:500]<-(-1) ; updown[(length(updown)-(500-1)):length(updown)]<-1 ; names(updown)<-names(tvalues) ; onesandzerosm[[nome[a]]]<-updown
}
genes<-names(onesandzerosm[[1]])
newmatrix<-c()
for(a in 1:length(names(onesandzerosm))){newmatrix<-rbind(newmatrix,onesandzerosm[[a]][genes]) ; print(paste(a," out of ",length(names(onesandzerosm)),sep=""))}
rownames(newmatrix)<-nome
save(newmatrix,file="Data/Matrix_with_patients_and_genes_1_0_-1.Rdata")
finishing<-Sys.time() ; print(finishing-starting)

print("End of part 1")

## Lets write the number of overlapping genes in each of the conditions! ##
## @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @ @ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ ##
load("Data/Top_500_sDEGs_names.Rdata")
fisin<-Sys.time()
for (a in 1:(length(names(degs))-1)){
  positi1<-degs[[a]]$up ; negati1<-degs[[a]]$down ; pega<-c() ; pegados<-c()
  for (b in (a+1):length(names(degs))){
    positi2<-degs[[b]]$up ; negati2<-degs[[b]]$down
    interspd <- length(intersect(positi1,positi2))
    intersnd <- length(intersect(negati1,negati2))
    interspi <- length(intersect(positi1,negati2))
    intersni <- length(intersect(negati1,positi2))
    junt<-c(names(degs)[a],names(degs)[b],interspi,interspd,intersni,intersnd)
    jun<-t(junt)
    write.table(jun,"Fisher_genes.txt",quote=F,sep="\t",row.names=F,col.names=F,append = T)
  }
  print(a/(length(names(degs))-1))
}
print("Finished!!") ; fisfin<-Sys.time() ; print(fisfin-fisin)



                                   ## @ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @ ##
                                   ## Calculate the number of intra-disease interactions increasing the number of selected genes ##
                                   ## @ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @ ##

if("Supplementary_documents"%in%list.files()==FALSE){dir.create("Supplementary_documents")}
if("Different_number_of_genes"%in%list.files("./Supplementary_documents")==FALSE){dir.create("Supplementary_documents/Different_number_of_genes")}
if("Plots"%in%list.files()==FALSE){dir.create("Plots")}

## Run lines from 334 to 387 assigning the following numbers to the object "selec" --> 100, 200, 300, 400, 500, 1000, 2000, 3000, 4000, 5000
starting<-Sys.time()
patients<-list.files("Data/Patienttvalues/") ; nome<-gsub(".Rdata","",patients,fixed=T) ; degs<-list()
## we select the number of genes as differentially expressed (100,200,300,400,500,1000,2000,3000,4000,5000)
number<-500
for(a in 1:length(patients)){
  load(paste("Data/Patienttvalues/",patients[a],sep=""))
  tvalues<-tvalues[order(as.numeric(tvalues),decreasing=T)]
  up<-names(tvalues)[1:number] ; down<-names(tvalues)[(length(tvalues)-(number-1)):length(tvalues)]
  degs[[gsub(".Rdata","",patients[a])]]$up<-up ; degs[[gsub(".Rdata","",patients[a])]]$down<-down
}
## Fisher test ##
background<-20502
diseases<-gsub("_.+","",names(degs)) ; disease<-unique(diseases)
for(z in 1:length(disease)){
  cojo<-which(diseases==disease[z])
  digs<-list() ; for(y in cojo){digs[[names(degs)[y]]]<-degs[y]}
  for (a in 1:(length(names(digs))-1)){
    positi1<-digs[[a]][[1]]$up ; negati1<-digs[[a]][[1]]$down ; pega<-c() ; pegados<-c()
    for (b in (a+1):length(names(digs))){
      positi2<-digs[[b]][[1]]$up ; negati2<-digs[[b]][[1]]$down
      interspd <- length(intersect(positi1,positi2))
      kkpd <- matrix(c(interspd,length(positi1)-interspd,length(positi2)-interspd,background+interspd-length(positi1)-length(positi2)),nrow=2,ncol=2)
      fispd<-fisher.test(kkpd,alternative="greater") $p.value
      intersnd <- length(intersect(negati1,negati2))
      kknd <- matrix(c(intersnd,length(negati1)-intersnd,length(negati2)-intersnd,background+intersnd-length(negati1)-length(negati2)),nrow=2,ncol=2)
      fisnd<-fisher.test(kknd,alternative="greater") $p.value
      interspi <- length(intersect(positi1,negati2))
      kkpi <- matrix(c(interspi,length(positi1)-interspi,length(negati2)-interspi,background+interspi-length(positi1)-length(negati2)),nrow=2,ncol=2)
      fispi<-fisher.test(kkpi,alternative="greater") $p.value
      intersni <- length(intersect(negati1,positi2))
      kkni <- matrix(c(intersni,length(negati1)-intersni,length(positi2)-intersni,background+intersni-length(negati1)-length(positi2)),nrow=2,ncol=2)
      fisni<-fisher.test(kkni,alternative="greater") $p.value
      junt<-c(names(digs)[a],names(digs)[b],fispi,fispd,fisni,fisnd)
      jun<-t(junt)
      write.table(jun,paste("Supplementary_documents/Different_number_of_genes/Fisher_overlaps_",number,".txt",sep=""),quote=F,sep="\t",row.names=F,col.names=F,append = T)
    }
    print(a/(length(names(digs))-1))
  }
}
ficheros<-list.files("Supplementary_documents/Different_number_of_genes/")
suma<-c()
for(a in 1:length(ficheros)){
  ## read the file and generate the network
  fis<-read.csv(paste("Supplementary_documents/Different_number_of_genes/",ficheros[a],sep=""),stringsAsFactors = F,sep="\t",header=F) ; threshold<-0.05
  negs<-list("Pos_Neg"=which(fis[,3]<=threshold),"Neg_Pos"=which(fis[,5]<=threshold),"Pos_Pos"=which(fis[,4]>threshold),"Neg_Neg"=which(fis[,6]>threshold))
  poss<-list("Pos_Neg"=which(fis[,3]>threshold),"Neg_Pos"=which(fis[,5]>threshold),"Pos_Pos"=which(fis[,4]<=threshold),"Neg_Neg"=which(fis[,6]<=threshold))
  interactions<-rep(0,length(fis[,1])) ; interactions[Reduce(intersect,negs)]<--1 ; interactions[Reduce(intersect,poss)]<-1
  tabint<-table(interactions) ; neg<-as.numeric(tabint[which(names(tabint)=="-1")]) ; pos<-as.numeric(tabint[which(names(tabint)=="1")]) ; suma<-c(suma,neg+pos)
}
seleccion<-as.numeric(gsub(".txt","",gsub("Fisher_overlaps_","",ficheros)))
pdf(file="Plots/Intra-disease_interactions_varying_number_of_deg.pdf")
plot(seleccion,suma,xlab="Number of Differentially Expressed Genes selected",ylab="Number of intra-disease interactions")
dev.off()


                                                ## @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@  @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ ##
                                                ## Extract the optimal threshold to obtain 0 interactions randomly ##
                                                ## @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@  @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ ##

## extract the analyzed genes
if ("Random_networks"%in%list.files("Data/") == FALSE){dir.create("Data/Random_networks")}
if ("Networks_varying_thresholds"%in%list.files("Data/") == FALSE){dir.create("Data/Networks_varying_thresholds")}
if ("Patient_similarities"%in%list.files("Data/Networks_varying_thresholds/") == FALSE){dir.create("Data/Networks_varying_thresholds/Patient_similarities")}
if ("Identifiers"%in%list.files("Data/Networks_varying_thresholds/") == FALSE){dir.create("Data/Networks_varying_thresholds/Identifiers")}
if ("RR"%in%list.files("Data/Networks_varying_thresholds/") == FALSE){dir.create("Data/Networks_varying_thresholds/RR")}
if ("Plots"%in%list.files("Data/Networks_varying_thresholds/") == FALSE){dir.create("Data/Networks_varying_thresholds/Plots")}
if ("Results"%in%list.files("Data/Networks_varying_thresholds/") == FALSE){dir.create("Data/Networks_varying_thresholds/Results")}
if ("Results"%in%list.files() == FALSE){dir.create("Results")}

## generate networks assigning randomly differentially expressed genes to each patient 
load("Data/Gene_Symbols_u133plus2.Rdata")
genes<-unique(simbolos)
set.seed(1)

## Run lines from 410 to 436 assigning the following numbers to the object "selec" --> 100, 200, 300, 400, 500, 1000, 2000, 3000, 4000, 5000
selec<-500
if(paste("Selecting_",selec,sep="")%in%list.files("Data/")==FALSE){dir.create(paste("Data/Random_networks_different_number_of_genes/",paste("Selecting_",selec,sep=""),sep=""))}
for(vuelta in 1:100){
  ldegs<-list() ; for(a in 1:1000){up<-sample(genes,selec) ; ldegs[[paste("Patient",a,sep="_")]]$up<-up ; ldegs[[paste("Patient",a,sep="_")]]$down<-sample(setdiff(genes,up),selec)}
  for (a in 1:(length(names(ldegs))-1)){
    positi1<-ldegs[[a]]$up ; negati1<-ldegs[[a]]$down ; pega<-c() ; pegados<-c()
    for (b in (a+1):length(names(ldegs))){
      positi2<-ldegs[[b]]$up ; negati2<-ldegs[[b]]$down
      interspd <- length(intersect(positi1,positi2))
      kkpd <- matrix(c(interspd,length(positi1)-interspd,length(positi2)-interspd,length(genes)+interspd-length(positi1)-length(positi2)),nrow=2,ncol=2)
      fispd<-fisher.test(kkpd,alternative="greater") $p.value
      intersnd <- length(intersect(negati1,negati2))
      kknd <- matrix(c(intersnd,length(negati1)-intersnd,length(negati2)-intersnd,length(genes)+intersnd-length(negati1)-length(negati2)),nrow=2,ncol=2)
      fisnd<-fisher.test(kknd,alternative="greater") $p.value
      interspi <- length(intersect(positi1,negati2))
      kkpi <- matrix(c(interspi,length(positi1)-interspi,length(negati2)-interspi,length(genes)+interspi-length(positi1)-length(negati2)),nrow=2,ncol=2)
      fispi<-fisher.test(kkpi,alternative="greater") $p.value
      intersni <- length(intersect(negati1,positi2))
      kkni <- matrix(c(intersni,length(negati1)-intersni,length(positi2)-intersni,length(genes)+intersni-length(negati1)-length(positi2)),nrow=2,ncol=2)
      fisni<-fisher.test(kkni,alternative="greater") $p.value
      junt<-c(names(ldegs)[a],names(ldegs)[b],fispi,fispd,fisni,fisnd)
      jun<-t(junt)
      write.table(jun,paste("Data/Random_networks_different_number_of_genes/",paste("Selecting_",selec,sep=""),"/Fisher_random_overlaps_",vuelta,".txt",sep=""),quote=F,sep="\t",row.names=F,col.names=F,append = T)
    }
  }
  print(paste("Round ",vuelta," finished in ",Sys.time()-inicio,sep=""))
}

## extract the number of random-interactions varying the fisher-test threshold
## Run lines from 440 to 463 assigning the following numbers to the object "selec" --> 100, 200, 300, 400, 500, 1000, 2000, 3000, 4000, 5000
selec<-500
thresholds<-c(0.05,0.01,0.005,0.001,0.0005,0.0001,0.00001) ; thressults<-list() ; means<-c() # we don't calculate 0.00005 because it involves the same number of genes (28)
for(threshold in thresholds){
  thressult<-c()
  for(a in 1:length(list.files(paste("Data/Random_networks_different_number_of_genes/Selecting_",selec,sep="")))){
    fis<-read.csv(paste(paste("Data/Random_networks_different_number_of_genes/Selecting_",selec,"/",sep=""),list.files(paste("Data/Random_networks_different_number_of_genes/Selecting_",selec,sep=""))[a],sep=""),stringsAsFactors = F,sep="\t",header=F)
    negs<-list("Pos_Neg"=which(fis[,3]<=threshold),"Neg_Pos"=which(fis[,5]<=threshold),"Pos_Pos"=which(fis[,4]>threshold),"Neg_Neg"=which(fis[,6]>threshold))
    poss<-list("Pos_Neg"=which(fis[,3]>threshold),"Neg_Pos"=which(fis[,5]>threshold),"Pos_Pos"=which(fis[,4]<=threshold),"Neg_Neg"=which(fis[,6]<=threshold))
    interactions<-rep(0,length(fis[,1])) ; interactions[Reduce(intersect,negs)]<--1 ; interactions[Reduce(intersect,poss)]<-1
    if(length(names(table(interactions)))>1){thressult<-c(thressult,sum(table(interactions)[which(names(table(interactions))!=0)]))}
    if(length(names(table(interactions)))==1){thressult<-c(thressult,0)} ; print(a)
  }
  thressults[[paste("Thres",threshold,sep="_")]]<-thressult
  means<-c(means,mean(thressult))
}
## represent a barplot showing the mean number of interactions
color<-rep("light blue",length(means))
color[which(means>0)]<-"red"
xx<-barplot(means,names.arg=thresholds,cex.names=0.8,cex.axis=0.8,ylim=c(0,max(thressults[[1]])+20),
            main="Mean number of random interactions\nvarying thresholds",col=color,border=color) ; text(x=xx,y=means,label=round(means,2),cex=0.8,pos=3,col=color)
pdf(paste("Plots/Mean_random_interactions_selecting_",selec,"_genes.pdf"))
xx<-barplot(means,names.arg=thresholds,cex.names=0.8,cex.axis=0.8,ylim=c(0,max(thressults[[1]])+20),
            main="Mean number of random interactions\nvarying thresholds",col=color,border=color) ; text(x=xx,y=means,label=round(means,2),cex=0.8,pos=3,col=color)
dev.off()


                                         ## @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @ @ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ ##
                                         ## To run this part you need to run first script 1 (to calculate fisher tests) ##
                                         ## @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @ @ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ ##
                                                              ## Generate ICD9, ICD10, and disease ##
                                                              ## @@ @@ @@ @@ @@ @ @ @@ @@ @@ @@ @@ ##
                                            ## Calculate patient-patient similarities based on different thresholds ##

## assign only 500 to the object "selec" at this point
selec<-500
if(paste("Selecting_",selec,"_genes",sep="")%in%list.files("Data/Networks_varying_thresholds/")==FALSE){dir.create(paste("Data/Networks_varying_thresholds/Selecting_",selec,"_genes",sep=""))}
if("Identifiers"%in%list.files(paste("Data/Networks_varying_thresholds/Selecting_",selec,"_genes",sep=""))==FALSE){dir.create(paste("Data/Networks_varying_thresholds/Selecting_",selec,"_genes/Identifiers",sep=""))}
if("Patient_similarities"%in%list.files(paste("Data/Networks_varying_thresholds/Selecting_",selec,"_genes",sep=""))==FALSE){dir.create(paste("Data/Networks_varying_thresholds/Selecting_",selec,"_genes/Patient_similarities",sep=""))}
if("Results"%in%list.files(paste("Data/Networks_varying_thresholds/Selecting_",selec,"_genes",sep=""))==FALSE){dir.create(paste("Data/Networks_varying_thresholds/Selecting_",selec,"_genes/Results",sep=""))}
if("RR"%in%list.files(paste("Data/Networks_varying_thresholds/Selecting_",selec,"_genes",sep=""))==FALSE){dir.create(paste("Data/Networks_varying_thresholds/Selecting_",selec,"_genes/RR",sep=""))}

fis<-read.csv(paste("Supplementary_documents/Different_number_of_genes/Fisher_overlaps_",selec,".txt",sep=""),stringsAsFactors = F,sep="\t",header=F)

# thresholds<-c(0.05,0.01,0.005,0.001,0.0005,0.0001,0.00001) # we don't calculate 0.00005 because it involves the same number of genes (28)
for(threshold in thresholds){
  negs<-list("Pos_Neg"=which(fis[,3]<=threshold),"Neg_Pos"=which(fis[,5]<=threshold),"Pos_Pos"=which(fis[,4]>threshold),"Neg_Neg"=which(fis[,6]>threshold))
  poss<-list("Pos_Neg"=which(fis[,3]>threshold),"Neg_Pos"=which(fis[,5]>threshold),"Pos_Pos"=which(fis[,4]<=threshold),"Neg_Neg"=which(fis[,6]<=threshold))
  interactions<-rep(0,length(fis[,1])) ; interactions[Reduce(intersect,negs)]<--1 ; interactions[Reduce(intersect,poss)]<-1
  ## save the patient-patient similarity network including 0s
  fisher_table<-cbind(fis[,1:2],interactions) ; colnames(fisher_table)<-c("Disease_1","Disease_2","Interaction")
  save(fisher_table,file=paste("Data/Networks_varying_thresholds/Selecting_",selec,"_genes/Patient_similarities/Patient_patient_similarity_network_",gsub(".","-",threshold,fixed=T),".Rdata",sep=""))
  print(which(thresholds==threshold))
}

## Generate the networks ##
## @@ @@ @@ @ @ @@ @@ @@ ##
# thresholds<-c(0.05,0.01,0.005,0.001,0.0005,0.0001,0.00001) # we don't calculate 0.00005 because it involves the same number of genes (28)
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
  sal<-c(length(which(cunpos>=length(intersect(posipeg,netwpeg))))/length(cunpos),length(which(cuneg>=length(intersect(negapeg,netwpeg))))/length(cuneg),
         length(intersect(posipeg,netwpeg)),length(intersect(negapeg,netwpeg)),length(posipeg),length(negapeg),length(netwpeg))
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
  }
  dev.off()
  results<-list("Epidemiological"=netwpeg,"Positive"=posipeg,"Negative"=negapeg,"Notconnected"=notconnected,"TotalNumbers"=totalnumbers,"Abstract"=sal)
  return(results)
}
ppsn<-list.files("Data/Networks_varying_thresholds/Patient_similarities/")[c(6,5,7,1,2,3,4)]

## ICD10 ##
## @@ @@ ##
load(paste("Data/Networks_varying_thresholds/Selecting_",selec,"_genes/Patient_similarities/",ppsn[1],sep="")) # fisher_table
patinf<-read.csv2("Patients_information_table.txt",stringsAsFactors = F,sep="\t",header=T)
splited<-strsplit(patinf[,1],"_")  ;  newnames<-c()
for(b in 1:length(splited)){newnames<-c(newnames,paste(patinf[b,11],splited[[b]][2],splited[[b]][3],splited[[b]][4],sep="_"))}
names(newnames)<-patinf[,1]
net<-cbind(as.character(newnames[as.character(fisher_table[,1])]),as.character(newnames[as.character(fisher_table[,2])]),as.character(fisher_table[,3]))
colnames(net)<-colnames(fisher_table) ; net0001<-net
disease<-unique(gsub("_.+","",unique(c(net0001[,1],net0001[,2]))))
one<-gsub("_.+","",net0001[,1],fixed=F)  ;  two<-gsub("_.+","",net0001[,2],fixed=F)  ;  intra<-which(one==two)  ;  inter0001<-net0001[-intra,]
## extract subgroup-identifiers
identificadores<-list()  ;  print("Starting the loop")
for(b in 1:length(disease)){
  ## select the identifiers of interest
  ini<-Sys.time();identificadores[[disease[b]]]<-c(grep(paste(disease[b],"_",sep=""),inter0001[,1]),grep(paste(disease[b],"_",sep=""),inter0001[,2]));fin<-Sys.time()
  print(paste(b," out of ",length(disease)," in ",fin-ini,sep=""))
}
save(identificadores,file=paste("Data/Networks_varying_thresholds/Selecting_",selec,"_genes/Identifiers/ICD10.Rdata",sep=""))
## calculate relative risks
for(a in 1:length(ppsn)){
  beggining<-Sys.time()
  load(paste("Data/Networks_varying_thresholds/Selecting_",selec,"_genes/Patient_similarities/",ppsn[a],sep="")) # fisher_table
  ## do the same as in the last step, just changing disease name by the corresponding ICD10 code
  patinf<-read.csv2("Patients_information_table.txt",stringsAsFactors = F,sep="\t",header=T)
  splited<-strsplit(patinf[,1],"_")  ;  newnames<-c()
  for(b in 1:length(splited)){newnames<-c(newnames,paste(patinf[b,11],splited[[b]][2],splited[[b]][3],splited[[b]][4],sep="_"))}
  names(newnames)<-patinf[,1]
  net<-cbind(as.character(newnames[as.character(fisher_table[,1])]),as.character(newnames[as.character(fisher_table[,2])]),as.character(fisher_table[,3]))
  colnames(net)<-colnames(fisher_table) ; net0001<-net
  disease<-unique(gsub("_.+","",unique(c(net0001[,1],net0001[,2]))))
  one<-gsub("_.+","",net0001[,1],fixed=F)  ;  two<-gsub("_.+","",net0001[,2],fixed=F)  ;  intra<-which(one==two)  ;  inter0001<-net0001[-intra,]
  ## Calculate Relative Risks ##
  load(paste("Data/Networks_varying_thresholds/Selecting_",selec,"_genes/Identifiers/ICD10.Rdata",sep=""))
  path<-paste("Data/Networks_varying_thresholds/Selecting_",selec,"_genes/RR/ICD10_",gsub(".Rdata","",strsplit(ppsn[a],"_")[[1]][5],fixed=T),sep="")
  estimate_rr(disease,identificadores,inter0001,path,1.96)
  info<-read.csv2("Patients_information_table.txt",stringsAsFactors = F,sep="\t",header=T)
  asoc<-info[,c(3,4,11)][-which(duplicated(info[,c(4,11)])),]
  ## Disease Trajectories - DMSN ##
  ## @@ @@ @@ @@ @ @ @@ @@ @@ @@ ##
  # Load our ICD10 interaction network and the disease trajectories, selecting only the significant interactions
  load(paste("Data/Networks_varying_thresholds/Selecting_",selec,"_genes/RR/ICD10_",gsub(".Rdata","",strsplit(ppsn[a],"_")[[1]][5],fixed=T),"_pRR.Rdata",sep=""))
  load(paste("Data/Networks_varying_thresholds/Selecting_",selec,"_genes/RR/ICD10_",gsub(".Rdata","",strsplit(ppsn[a],"_")[[1]][5],fixed=T),"_nRR.Rdata",sep=""))
  
  network<-read.csv2("Epidemiological_networks/Disease_pairs.csv",stringsAsFactors = F,sep=",",header=T)[,1:3]
  dtp<-network_overlap(RRpos,RRneg,network,
                       paste("Data/Networks_varying_thresholds/Selecting_",selec,"_genes/Plots/Pairs_DTP-DMSN_",gsub(".Rdata","",strsplit(ppsn[a],"_")[[1]][5],fixed=T),".pdf",sep=""))
  save(dtp,file=paste("Data/Networks_varying_thresholds/Selecting_",selec,"_genes/Results/Summary_ICD10_",gsub(".Rdata","",strsplit(ppsn[a],"_")[[1]][5],fixed=T),".Rdata",sep=""))
  print(paste("Round ",a," finished",sep=""))
  print(Sys.time()-beggining)
}

## ICD9 ##
## @@@@ ##
load(paste("Data/Networks_varying_thresholds/Selecting_",selec,"_genes/Patient_similarities/",ppsn[1],sep="")) # fisher_table
## do the same as in the last step, just changing disease name by the corresponding ICD10 code
patinf<-read.csv2("Patients_information_table.txt",stringsAsFactors = F,sep="\t",header=T)
splited<-strsplit(patinf[,1],"_")  ;  newnames<-c()
for(b in 1:length(splited)){newnames<-c(newnames,paste(gsub("_","-",patinf[b,4]),splited[[b]][2],splited[[b]][3],splited[[b]][4],sep="_"))}
names(newnames)<-patinf[,1]
net<-cbind(as.character(newnames[as.character(fisher_table[,1])]),as.character(newnames[as.character(fisher_table[,2])]),as.character(fisher_table[,3]))
colnames(net)<-colnames(fisher_table) ; net0001<-net ; disease<-unique(gsub("_.+","",unique(c(net0001[,1],net0001[,2]))))
one<-gsub("_.+","",net0001[,1],fixed=F)  ;  two<-gsub("_.+","",net0001[,2],fixed=F)  ;  intra<-which(one==two)  ;  inter0001<-net0001[-intra,]
## extract subgroup-identifiers
identificadores<-list()  ;  print("Starting the loop")
for(b in 1:length(disease)){
  ini<-Sys.time();identificadores[[disease[b]]]<-c(grep(paste(disease[b],"_",sep=""),inter0001[,1]),grep(paste(disease[b],"_",sep=""),inter0001[,2]));fin<-Sys.time()
  print(paste(b," out of ",length(disease)," in ",fin-ini,sep=""))
}
save(identificadores,file=paste("Data/Networks_varying_thresholds/Selecting_",selec,"_genes/Identifiers/ICD9.Rdata",sep=""))
## calculate relative risks
for(a in 1:length(ppsn)){
  beggining<-Sys.time()
  load(paste("Data/Networks_varying_thresholds/Selecting_",selec,"_genes/Patient_similarities/",ppsn[a],sep="")) # fisher_table
  ## do the same as in the last step, just changing disease name by the corresponding ICD10 code
  patinf<-read.csv2("Patients_information_table.txt",stringsAsFactors = F,sep="\t",header=T)
  splited<-strsplit(patinf[,1],"_")  ;  newnames<-c()
  for(b in 1:length(splited)){newnames<-c(newnames,paste(gsub("_","-",patinf[b,4]),splited[[b]][2],splited[[b]][3],splited[[b]][4],sep="_"))}
  names(newnames)<-patinf[,1]
  net<-cbind(as.character(newnames[as.character(fisher_table[,1])]),as.character(newnames[as.character(fisher_table[,2])]),as.character(fisher_table[,3]))
  colnames(net)<-colnames(fisher_table) ; net0001<-net
  disease<-unique(gsub("_.+","",unique(c(net0001[,1],net0001[,2]))))
  one<-gsub("_.+","",net0001[,1],fixed=F)  ;  two<-gsub("_.+","",net0001[,2],fixed=F)  ;  intra<-which(one==two)  ;  inter0001<-net0001[-intra,]
  load(paste("Data/Networks_varying_thresholds/Selecting_",selec,"_genes/Identifiers/ICD9.Rdata",sep=""))
  ## Calculate Relative Risks ##
  path<-paste("./Data/Networks_varying_thresholds/Selecting_",selec,"_genes/RR/ICD9_",gsub(".Rdata","",strsplit(ppsn[a],"_")[[1]][5],fixed=T),sep="")
  ## since this network is being generated to compare it with the PDN we use a more restrictive threshold (the one used by Hidalgo et al. 2009)
  estimate_rr(disease,identificadores,inter0001,path,2.576)
  info<-read.csv2("Patients_information_table.txt",stringsAsFactors = F,sep="\t",header=T)
  ## Disease Trajectories - DMSN ##
  ## @@ @@ @@ @@ @ @ @@ @@ @@ @@ ##
  # Load our ICD10 interaction network and the disease trajectories, selecting only the significant interactions
  load(paste("Data/Networks_varying_thresholds/Selecting_",selec,"_genes/RR/ICD9_",gsub(".Rdata","",strsplit(ppsn[a],"_")[[1]][5],fixed=T),"_pRR.Rdata",sep=""))
  load(paste("Data/Networks_varying_thresholds/Selecting_",selec,"_genes/RR/ICD9_",gsub(".Rdata","",strsplit(ppsn[a],"_")[[1]][5],fixed=T),"_nRR.Rdata",sep=""))
  husig<-read.csv2("Epidemiological_networks/PDN_3_digits.net",stringsAsFactors = F,header=F,sep="\t")
  network2<-cbind(paste("ICD9-",husig[,1],sep=""),paste("ICD9-",husig[,2],sep=""),husig[,7])
  network<-rbind(cbind(as.character(network2[,1]),as.character(network2[,2]),network2[,3]),
                 cbind(as.character(network2[,2]),as.character(network2[,1]),network2[,3]))
  dtp<-network_overlap(RRpos,RRneg,network,
                       paste("Data/Networks_varying_thresholds/Selecting_",selec,"_genes/Plots/Pairs_PDN-DMSN_",gsub(".Rdata","",strsplit(ppsn[a],"_")[[1]][5],fixed=T),".pdf",sep=""))
  save(dtp,file=paste("Data/Networks_varying_thresholds/Selecting_",selec,"_genes/Results/Summary_ICD9_",gsub(".Rdata","",strsplit(ppsn[a],"_")[[1]][5],fixed=T),".Rdata",sep=""))
  print(paste("Round ",a," finished",sep=""))
  print(Sys.time()-beggining)
}

## Disease ##
## @@ @ @@ ##
load(paste("Data/Networks_varying_thresholds/Selecting_",selec,"_genes/Patient_similarities/",ppsn[1],sep="")) ; net0001<-fisher_table
disease<-unique(gsub("_.+","",unique(c(net0001[,1],net0001[,2]))))
one<-gsub("_.+","",net0001[,1],fixed=F)  ;  two<-gsub("_.+","",net0001[,2],fixed=F)  ;  intra<-which(one==two)  ;  inter0001<-net0001[-intra,]
identificadores<-list()  ;  print("Starting the loop")
for(b in 1:length(disease)){
  ini<-Sys.time();identificadores[[disease[b]]]<-c(grep(paste(disease[b],"_",sep=""),inter0001[,1]),grep(paste(disease[b],"_",sep=""),inter0001[,2]));fin<-Sys.time()
  print(paste(b," out of ",length(disease)," in ",fin-ini,sep=""))
}
save(identificadores,file=paste("Data/Networks_varying_thresholds/Selecting_",selec,"_genes/Identifiers/Disease.Rdata",sep=""))
for(a in 1:length(ppsn)){
  beggining<-Sys.time()
  load(paste("Data/Networks_varying_thresholds/Selecting_",selec,"_genes/Patient_similarities/",ppsn[a],sep="")) # fisher_table
  net0001<-fisher_table
  disease<-unique(gsub("_.+","",unique(c(net0001[,1],net0001[,2]))))
  one<-gsub("_.+","",net0001[,1],fixed=F)  ;  two<-gsub("_.+","",net0001[,2],fixed=F)  ;  intra<-which(one==two)  ;  inter0001<-net0001[-intra,]
  load(paste("Data/Networks_varying_thresholds/Selecting_",selec,"_genes/Identifiers/Disease.Rdata",sep=""))
  ## Calculate Relative Risks ##
  path<-paste("./Data/Networks_varying_thresholds/Selecting_",selec,"_genes/RR/Disease_",gsub(".Rdata","",strsplit(ppsn[a],"_")[[1]][5],fixed=T),sep="")
  estimate_rr(disease,identificadores,inter0001,path,1.96)
  print(paste("Round ",a," finished",sep="")) ; print(Sys.time()-beggining)
}

## Load the overlap results and generate a summary-table ##
## @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ ##
## ICD9 ##
## @@@@ ##
results<-list.files(paste("Data/Networks_varying_thresholds/Selecting_",selec,"_genes/Results/",sep="")) ; results<-results[grep("ICD9",results)] ; traject_overlap9<-c()
for(a in 1:length(results)){
  load(paste("Data/Networks_varying_thresholds/Selecting_",selec,"_genes/Results/",results[a],sep="")) # dtp
  traject_overlap9<-rbind(traject_overlap9,c(gsub(".Rdata","",gsub("Summary_ICD9_","",results[a]),fixed=T),dtp$Abstract))
}
colnames(traject_overlap9)<-c("Threshold","9+pval","9-pval","9i+n","9i-n","9+n","9-n","9tn") ; traject_overlap9<-traject_overlap9[c(4,3,2,1,7,5,6),]
## ICD10 ##
## @@ @@ ##
results<-list.files(paste("Data/Networks_varying_thresholds/Selecting_",selec,"_genes/Results/",sep="")) ; results<-results[grep("ICD10",results)] ; traject_overlap10<-c()
for(a in 1:length(results)){
  load(paste("Data/Networks_varying_thresholds/Selecting_",selec,"_genes/Results/",results[a],sep="")) # dtp
  traject_overlap10<-rbind(traject_overlap10,c(gsub(".Rdata","",gsub("Summary_ICD10_","",results[a]),fixed=T),dtp$Abstract))
}
colnames(traject_overlap10)<-c("Threshold","10+pval","10-pval","10i+n","10i-n","10+n","10-n","10tn") ; traject_overlap10<-traject_overlap10[c(4,3,2,1,7,5,6),]
traject_overlap<-cbind(traject_overlap9,traject_overlap10[,2:8])
write.table(traject_overlap,paste("Results/Significance_of_the_overlaps_varying_thresholds_selecting_",selec,"_genes.txt",sep=""),quote=F,sep="\t",row.names=F)


