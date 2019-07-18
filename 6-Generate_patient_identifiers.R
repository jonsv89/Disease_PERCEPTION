## @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @ @ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ ##
## We generate the identifiers of each patient to generate the patient-specific comorbidity profile ##
## @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @ @ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ ##

if("Patient_identifiers"%in%list.files("Data/Identifiers")==FALSE){dir.create("Data/Identifiers/Patient_identifiers")}
load("Data/Patient_patient_similarity_network.Rdata")  ;  red000005<-fisher_table

enfermedad<-unique(gsub("_.+","",unique(c(red000005[,1],red000005[,2]))))
## Remove intra-disease interactions
uno<-gsub("_.+","",red000005[,1],fixed=F)  ;  dos<-gsub("_.+","",red000005[,2],fixed=F)  ;  intra<-which(uno==dos)  ;  inter000005<-red000005[-intra,]
tablon<-read.csv2("Patients_information_table.txt",stringsAsFactors = F,sep="\t",header=T)
pacientes<-tablon[,1]
## These are the arguments:
seq(1,136,28)
args=commandArgs(trailingOnly = TRUE)  ;  argumento<-as.numeric(args)  ;  print(class(argumento))  ;  if(argumento==113){fin<-23}    if(argumento!=113){fin<-27}
## Start extracting the identifiers
for(b in argumento:(argumento+fin)){
  pacientos<-pacientes[grep(enfermedad[b],pacientes)]  ;  identificadores<-list()
  for(a in 1:length(pacientos)){ini<-Sys.time()  ;  identificadores[[pacientos[a]]]<-c(grep(pacientos[a],inter000005[,1]),grep(pacientos[a],inter000005[,2]))  ;  fin<-Sys.time()}
  print(paste(b," de ",length(enfermedad)," ",fin-ini,sep=""))
  save(identificadores,file=paste("Data/Identifiers/Patient_identifiers/",enfermedad[b],".RData",sep=""))
}
print("Finished!")