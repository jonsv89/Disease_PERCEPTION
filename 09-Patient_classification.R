## All the results are going to be saved in "Gene_based_patient_classification" ##
## @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@@@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ ##

if("Gene_based_patient_classification"%in%list.files("Data/")==FALSE){dir.create("./Data/Gene_based_patient_classification/")}

## Load functions ##
## @@ @@ @@ @@ @@ ##

## Check wether the patients properly assigned to their diseases are properly classified into their subgroup
extract_subgroup_classification<-function(gcasdif,gcasificados,pati){
  elprimero<-c()
  elcompleto<-c()
  es<-c()
  sale<-c()
  elpat<-c()
  for(a in 1:length(names(gcasdif))){
    enfpat<-gsub("_.+","",names(gcasdif)[a])
    asigned_disease<-gcasdif[[a]][1]
    if(enfpat==asigned_disease){
      print(paste(names(gcasdif)[a],"OK!",sep="              "))
      print(a)
      grupos<-as.character(pati[names(gcasificados[[a]])])
      groups<-grupos[grep(enfpat,grupos)]
      elmio<-as.character(pati[names(gcasdif)[a]])
      # Fin the most similar patient #
      if(elmio==groups[1]){
        elprimero<-c(elprimero,"si")
        es<-c(es,elmio)
        sale<-c(sale,groups[1])
      }
      if(elmio!=groups[1]){
        elprimero<-c(elprimero,"no")
        es<-c(es,elmio)
        sale<-c(sale,groups[1])
      }
      # Check the first completed patient-subgroup
      grupis<-unique(groups)
      final<-c()
      for(b in 1:length(grupis)){
        final<-c(final,max(which(groups==grupis[b])))
      }
      names(final)<-grupis
      laclas<-names(which(final==min(final)))
      
      if(elmio==laclas){
        elcompleto<-c(elcompleto,"si")
      }
      if(elmio!=laclas){
        elcompleto<-c(elcompleto,"no")
      }
      elpat<-c(elpat,names(gcasdif)[a])
    }
    if(enfpat!=asigned_disease){
      elprimero<-c(elprimero,"-")
      elcompleto<-c(elcompleto,"-")
    }
  }
  resultado<-list("first"=elprimero,"complete"=elcompleto,"correct"=es,"result"=sale,"patient_name"=elpat)
  return(resultado)
}
## Sacar el numero de pacientes que clasifican en cada enfermedad
number_patients<-function(icasdif,num){
  cogio<-c()
  for(a in 1:length(names(icasdif))){
    seleccion<-names(which(table(icasdif[[a]][1:num])==max(table(icasdif[[a]][1:num]))))
    if(length(seleccion)>1){cogio<-c(cogio,"Varios")}
    if(length(seleccion)==1){cogio<-c(cogio,seleccion)}
  }
  return(cogio)
}
## Load necessary information ##
## @@ @@ @@ @@ @@ @@ @@ @@ @@ ##
## Patients information
pacientes<-read.csv2("Patients_information_table.txt",stringsAsFactors = F,sep="\t",header=T)
## Subgroups
pati<-pacientes$Cluster_info_cluster_specific
names(pati)<-pacientes$Patient
## Diseases
patid<-pacientes$Disease
names(patid)<-pacientes$Patient
## Colors
patic<-pacientes$Disease_group_color
names(patic)<-pacientes$Disease
                                                                ## ~ @@ @@ @@@@ @@ @@ ~ ##
                                                                ## ~ Lets get started ~ ##
                                                                ## ~ @@ @@ @@@@ @@ @@ ~ ##


## Check how patients are classified into their corresponding disease ##
## @@ @@ @@ @@ @@ @@ @@ @@ @ @ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ ##

load("Data/Matrix_with_patients_and_genes_1_0_-1.Rdata")
patgen<-t(newmatrix)
gcasificados<-list() ; gcasdif<-c()
nuestros<-patgen[,c(grep("Alzheimer",colnames(patgen)),grep("NSCLC",colnames(patgen)))]
for(a in 1:length(nuestros[1,])){
  paciente<-as.numeric(nuestros[,a]) ; diferencias<-c() ; patgen2<-patgen[,-which(colnames(patgen)==colnames(nuestros)[a])]
  for(b in 1:length(patgen2[1,])){compa<-as.numeric(patgen2[,b]) ; diferencias<-c(diferencias,dist(rbind(paciente,compa)))}
  names(diferencias)<-colnames(patgen2)
  diferencias<-sort(diferencias) ; gcasificados[[colnames(nuestros)[a]]]<-diferencias ; gcasdif[[colnames(nuestros)[a]]]<-gsub("_.+","",names(diferencias))
  print(gsub("_.+","",names(diferencias))[1])
}
save(gcasificados,gcasdif,file="Data/Gene_based_patient_classification/LOO_patient_clasification_euclidean_distance_AD_NSCLC.Rdata")
load("Data/Gene_based_patient_classification/LOO_patient_clasification_euclidean_distance_AD_NSCLC.Rdata")
## Results
tenim<-c();for(a in 1:length(names(gcasdif))){tenim<-c(tenim,gcasdif[[a]][1])}
as.data.frame(table(tenim[grep("Alzheimer",names(gcasdif))]))[order(as.data.frame(table(tenim[grep("Alzheimer",names(gcasdif))]))[,2],decreasing=T),]
as.data.frame(table(tenim[grep("NSCLC",names(gcasdif))]))[order(as.data.frame(table(tenim[grep("NSCLC",names(gcasdif))]))[,2],decreasing=T),]


## Check how patients are classified into their corresponding subgroup ##
## @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ ##

resultado<-extract_subgroup_classification(gcasdif,gcasificados,pati)
## Selecting the first one
table(resultado$first[grep("Alzheimer",names(gcasdif))])
table(resultado$first[grep("NSCLC",names(gcasdif))])


adno<-names(gcasdif)[grep("Alzheimer",names(gcasdif))][which(resultado$first[grep("Alzheimer",names(gcasdif))]=="no")]
adsi<-names(gcasdif)[grep("Alzheimer",names(gcasdif))][which(resultado$first[grep("Alzheimer",names(gcasdif))]=="si")]
lcno<-names(gcasdif)[grep("NSCLC",names(gcasdif))][which(resultado$first[grep("NSCLC",names(gcasdif))]=="no")]
lcsi<-names(gcasdif)[grep("NSCLC",names(gcasdif))][which(resultado$first[grep("NSCLC",names(gcasdif))]=="si")]

info<-read.csv2("Patients_information_table.txt",stringsAsFactors = F,sep="\t",header=T) ; clusclus<-info$Cluster_info_cluster_specific ; names(clusclus)<-info$Patient

## Percentage of patient not properly classified into their group associated to non-trustable subgroups (AD)
length(grep("NS$",as.character(clusclus[adno])))/length(as.character(clusclus[adno]))
length(grep("NS$",as.character(clusclus[adsi])))/length(as.character(clusclus[adsi]))
## Percentage of patient not properly classified into their group associated to non-trustable subgroups (NSCLC)
length(grep("NS$",as.character(clusclus[lcno])))/length(as.character(clusclus[lcno]))
length(grep("NS$",as.character(clusclus[lcsi])))/length(as.character(clusclus[lcsi]))


## Calculate for each subgroup the number of TP, FP,TN, FN ##
## @@ @@ @@ @@ @@ @@ @@ @@ @@ @ @@ @@ @@ @@ @@ @@ @@ @@ @@ ##

ecasdif<-list() ; for(a in 1:length(names(gcasificados))){ecasdif[[names(gcasificados)[a]]]<-as.character(pati[names(gcasificados[[a]])])}


disease<-"Alzheimer"
## Alzheimer case
subgrupos<-unique(as.character(pati[names(ecasdif)[(grep(disease,names(ecasdif)))]]))
tp<-c();fn<-c();tn<-c()
aparecen_cuando_no<-c()
for(a in 1:length(subgrupos)){
  salidas<-resultado$result[which(resultado$correct==subgrupos[a])]
  tp<-c(tp,length(which(salidas==subgrupos[a])))
  fn<-c(fn,length(which(salidas!=subgrupos[a])))
  aparecen_cuando_no<-rbind(aparecen_cuando_no,cbind(names(table(salidas[which(salidas!=subgrupos[a])])),as.numeric(table(salidas[which(salidas!=subgrupos[a])]))))
  result<-resultado$result[grep(disease,resultado$correct)]
  correct<-resultado$correct[grep(disease,resultado$correct)]
  tn<-c(tn,length(which(result[which(correct!=subgrupos[a])]!=subgrupos[a])))
}
names(tp)<-subgrupos
names(fn)<-subgrupos
names(tn)<-subgrupos
tienen<-unique(aparecen_cuando_no[,1])
fp<-c();
for(a in 1:length(tienen)){
  bu<-as.numeric(which(aparecen_cuando_no[,1]==tienen[a]))
  if(length(bu)==1){buf<-as.numeric(aparecen_cuando_no[bu,2]) ; names(buf)<-aparecen_cuando_no[bu,1] ; fp<-c(fp,buf)}
  if(length(bu)>1){buf<-sum(as.numeric(aparecen_cuando_no[bu,2])) ; names(buf)<-tienen[a] ; fp<-c(fp,buf)}
}
cero<-rep(0,length(setdiff(names(tp),names(fp))))
names(cero)<-setdiff(names(tp),names(fp))
fp<-c(fp,cero)
fp<-fp[subgrupos]
## Sensitivity
sensitivity<-as.numeric(tp)/(as.numeric(tp)+as.numeric(fn))
names(sensitivity)<-subgrupos
tamanos<-table(as.character(pati[names(ecasdif)[grep(disease,names(ecasdif))]]))[subgrupos]
plot(as.numeric(tamanos),as.numeric(sensitivity),xlab="Cluster size",ylab="Sensitivity",main=disease)
correlation<-cor.test(as.numeric(tamanos),as.numeric(sensitivity));correlation
abline(lm(as.numeric(sensitivity)~as.numeric(tamanos)),col="red")
text(10,y=0.7,paste("cor =",correlation$estimate,sep=" "))
## Specificity
specificity<-as.numeric(tn)/(as.numeric(tn)+as.numeric(fp))
names(sensitivity)<-subgrupos
tamanos<-table(as.character(pati[names(ecasdif)[grep(disease,names(ecasdif))]]))[subgrupos]
plot(as.numeric(tamanos),as.numeric(specificity),xlab="Cluster size",ylab="Specificity",main=disease,ylim=c(0.7,1))
correlation<-cor.test(as.numeric(tamanos),as.numeric(specificity));correlation
abline(lm(as.numeric(specificity)~as.numeric(tamanos)),col="red")
text(10,y=0.9,paste("cor =",correlation$estimate,sep=" "))
sensitivity/specificity
mcc<-(as.numeric(tp)*as.numeric(tn)-as.numeric(fp)*as.numeric(fn))/sqrt((as.numeric(tp)+as.numeric(fp))*(as.numeric(tp)+as.numeric(fn))*(as.numeric(tn)+as.numeric(fp))*(as.numeric(tn)+as.numeric(fn))) 

## @@ @@ @@ @@ @@ @@ @@ @@ @@ @@@@ @@ @@ @@ @@ @@ @@ @@ @@ @@ ##
## Classify all the patients and not only the AD - NSCLC ones ##
## @@ @@ @@ @@ @@ @@ @@ @@ @@ @@@@ @@ @@ @@ @@ @@ @@ @@ @@ @@ ##
gencasificados<-list()
gencasdif<-c()
load("Data/Matrix_with_patients_and_genes_1_0_-1.Rdata")
patgen<-t(newmatrix)
for(a in 1:length(patgen[1,])){
  paciente<-as.numeric(patgen[,a])
  diferencias<-c()
  patgen2<-patgen[,-a]
  for(b in 1:length(patgen2[1,])){compa<-as.numeric(patgen2[,b]) ; diferencias<-c(diferencias,dist(rbind(paciente,compa)))}
  names(diferencias)<-colnames(patgen2)
  diferencias<-sort(diferencias)
  gencasificados[[colnames(patgen)[a]]]<-diferencias
  gencasdif[[colnames(patgen)[a]]]<-gsub("_.+","",names(diferencias))
  print(gsub("_.+","",names(diferencias))[1])
  if(length(which(seq(0,6200,100)==a))>0){
    save(gencasificados,gencasdif,file="Data/Gene_based_patient_classification/LOO_patient_clasification_euclidean_distance_ALL.Rdata")
    print(paste("Hemos guardado",a,"casos",sep="  "))
  }
}
save(gencasificados,gencasdif,file="Data/Gene_based_patient_classification/LOO_patient_clasification_euclidean_distance_ALL.Rdata")
print("Terminado")

## Calculate false positives in our Alzheimer's disease analyses ##
## @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ ##
load("Data/Gene_based_patient_classification/LOO_patient_clasification_euclidean_distance_ALL.Rdata") # gencasificados gencasdif
## Alzheimer's disease ##
enfermedades<-unique(gsub("_.+","",names(gencasdif)))
cuantosporenfermedad<-as.data.frame(table(gsub("_.+","",names(gencasdif))))
# mccs<-c()
for(a in 1:length(enfermedades)){
  # a<-1
  primeraopcion<-number_patients(gencasdif,1)
  primeraopcion2<-primeraopcion[-grep(enfermedades[a],names(gencasdif))]
  fp<-length(which(primeraopcion2==enfermedades[a]))
  tn<-length(primeraopcion2)-fp
  cuantos<-primeraopcion[grep(enfermedades[a],names(gencasdif))]
  tp<-length(which(cuantos==enfermedades[a]))
  fn<-length(cuantos)-tp
  plantar<-c(tp,fp,tn,fn)
  write.table(t(plantar),"Data/Gene_based_patient_classification/TP_FP_TN_FN_values.txt",quote=F,row.names=F,col.names=F,append = T)
  ## Sensitivities and Specificities ##
  sensi<-c(tp/(tp+fn),tn/(tn+fp))
  write.table(t(sensi),"Data/Gene_based_patient_classification/Sensitivity_and_specificity_values.txt",quote=F,row.names=F,col.names=F,append = T)
  print(a)
}
sensi_speci<-read.table("Data/Gene_based_patient_classification/Sensitivity_and_specificity_values.txt",stringsAsFactors = F,sep=" ",header=F)
sensitivities<-sensi_speci[,1]
specificities<-sensi_speci[,2]


colores<-as.character(patic[enfermedades])
plot(sensitivities,specificities,col=colores)
rocno<-cbind(sensitivities,specificities)
rownames(rocno)<-enfermedades
colnames(rocno)<-c("Sensitivity","Specificity")
rocno<-rocno[order(rocno[,1],decreasing=T),]

## Percentage of hits ##
pred<-c() ; for(a in 1:length(names(gencasdif))){pred<-c(pred,gencasdif[[a]][1])}
which(gsub("_.+","",names(gencasdif))==pred)

## Plot for each disease the specificity vs. sensitivity ##
## @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ ##
values<-list()
sensitivities<-c()
specificities<-c()
sdspecificities<-c()
for(a in 1:length(enfermedades)){
  primeraopcion<-number_patients(gencasdif,1)
  primeraopcion2<-primeraopcion[-grep(enfermedades[a],names(gencasdif))]
  fp<-length(which(primeraopcion2==enfermedades[a]))
  cuantos<-primeraopcion[grep(enfermedades[a],names(gencasdif))]
  tp<-length(which(cuantos==enfermedades[a]))
  fn<-length(cuantos)-tp
  sensitivities<-c(sensitivities,tp/(tp+fn))
  if(fp==0){
    tn<-length(primeraopcion2)-fp
    specificities<-c(specificities,tn/(tn+fp))
    sdspecificities<-c(sdspecificities,0)
  }
  if(fp>0){
    fps<-c();tns<-c()
    for(b in 1:100){
      fip<-length(which(primeraopcion2[sample(1:length(primeraopcion2),length(cuantos))]==enfermedades[a]))
      fps<-c(fps,fip)
      tns<-c(tns,length(primeraopcion2[sample(1:length(primeraopcion2),length(cuantos))])-fip)
    }
    fp<-mean(fps)
    fpsd<-sd(fps)
    tn<-mean(tns)
    specificities<-c(specificities,mean(tns/(tns+fps)))
    sdspecificities<-c(sdspecificities,sd(tns/(tns+fps)))
  }
  values[[enfermedades[a]]]$tp<-tp
  values[[enfermedades[a]]]$fp<-fp
  values[[enfermedades[a]]]$tn<-tn
  values[[enfermedades[a]]]$fn<-fn
}
pdf(file="Plots/ROC_class_into_disease.pdf")
colores<-as.character(patic[enfermedades])
plot(sensitivities,specificities,col=colores,ylab="Specificity",xlab="Sensitivity")
dev.off()

rocno<-cbind(sensitivities,specificities)
rownames(rocno)<-enfermedades
colnames(rocno)<-c("Sensitivity","Specificity")
rocno<-rocno[order(rocno[,1],decreasing=T),]
tps<-c();for(a in 1:length(names(values))){tps<-c(tps,values[[a]]$tp)}
sum(tps)/length(names(gencasificados))

## Repeat the analisis generating pacients in a random way ##
## @ @ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @ @ ##
## Clasificar todos los pacientes en base a los genes ##
gencasificados<-list()
gencasdif<-c()
load("Data/Matrix_with_patients_and_genes_1_0_-1.Rdata")
patgen<-t(newmatrix)
for(a in 1:length(patgen[1,])){
  paciente<-as.numeric(patgen[,a])
  paciente<-sample(paciente)
  diferencias<-c()
  patgen2<-patgen[,-a]
  for(b in 1:length(patgen2[1,])){
    compa<-as.numeric(patgen2[,b])
    diferencias<-c(diferencias,dist(rbind(paciente,compa)))
  }
  names(diferencias)<-colnames(patgen2)
  diferencias<-sort(diferencias)
  gencasificados[[colnames(patgen)[a]]]<-diferencias
  gencasdif[[colnames(patgen)[a]]]<-gsub("_.+","",names(diferencias))
  print(gsub("_.+","",names(diferencias))[1])
  if(length(which(seq(0,6200,100)==a))>0){
    save(gencasificados,gencasdif,file="Data/Gene_based_patient_classification/LOO_patient_clasification_euclidean_distance_ALL_random.Rdata")
    print(paste("We have saved",a,"cases",sep="  "))
  }
}
save(gencasificados,gencasdif,file="Data/Gene_based_patient_classification/LOO_patient_clasification_euclidean_distance_ALL_random.Rdata")
print("Ended!")
