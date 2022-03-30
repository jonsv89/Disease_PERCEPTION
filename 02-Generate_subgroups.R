## @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @ @ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ ##
## Extract the optimal number of clusters (patient subgroups) per disease based on differential gene expression information ##
## @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @ @ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ ##
## estimated duration:
library("cluster");library("NbClust");library(e1071)
set.seed(1)
extract_subgroups<-function(matinter,enfermedad){
  if("Genes_based_clustering"%in%list.files("Data/")==FALSE){dir.create("Data/Genes_based_clustering")}
  if("Silhouette"%in%list.files("Data/Genes_based_clustering")==FALSE){dir.create("Data/Genes_based_clustering/Silhouette")}
  if("Clusters"%in%list.files("Data/Genes_based_clustering")==FALSE){dir.create("Data/Genes_based_clustering/Clusters")}
  k.max <- length(rownames(matinter))-1  ## Number of patients -1
  if(k.max>1){
    data <- matinter
    sil <- rep(0, k.max)
    # Compute the average silhouette width for
    # k = 2 to k = 15
    for(i in 2:k.max){
      tryCatch({km.res <-try( kmeans(data, centers = i, nstart = 100))
      ss <- silhouette(km.res$cluster, dist(data))
      sil[i] <- mean(ss[, 3])
      }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
    }
    # Plot the  average silhouette width
    plot(1:k.max, sil, type = "b", pch = 19, frame = FALSE, xlab = "Number of clusters k",main=enfermedad)
    abline(v = which.max(sil), lty = 2)
    pdf(file=paste("Data/Genes_based_clustering/Silhouette/",enfermedad,".pdf",sep=""))
    plot(1:k.max, sil, type = "b", pch = 19, frame = FALSE, xlab = "Number of clusters k",main=enfermedad)
    abline(v = which.max(sil), lty = 2)
    dev.off()
    km <- kmeans(matinter, centers = which.max(sil), nstart = 100)
    return(km)
  }
}

## We provide intervals (lines 49-53) for running the script from lines 36 to 67 in 5 cores at the same time
args=commandArgs(trailingOnly = TRUE)
argumento<-as.numeric(args)
print(class(argumento))

starting<-Sys.time()
load("Data/Matrix_with_patients_and_genes_1_0_-1.Rdata") #newmatrix
diseases<-gsub("_.+","",rownames(newmatrix))
disease<-unique(diseases)
info<-read.csv2("Patients_information_table.txt",stringsAsFactors = F,sep="\t",header=T)
numbers<-table(info[,3])
numbers<-numbers[order(numbers,decreasing = F)]
disease<-names(numbers)

if(argumento==1){sequences<-seq(1,135,5)}
if(argumento==2){sequences<-seq(2,135,5)}
if(argumento==3){sequences<-seq(3,135,5)}
if(argumento==4){sequences<-seq(4,135,5)}
if(argumento==5){sequences<-seq(5,135,5)}

print(sequences)
print("We start the loop!")
for(a in sequences){
  star<-Sys.time()
  dis<-disease[a]
  matinter<-newmatrix[grep(paste("^",dis,"$",sep=""),diseases),]
  ends<-Sys.time()
  clusters<-extract_subgroups(matinter,dis)
  save(clusters,file=paste("Data/Genes_based_clustering/Clusters/",dis,".Rdata",sep=""))
  print(paste(dis," finished, ",a," out of ",length(disease)," in ",ends-star,sep=" "))
}
finishing<-Sys.time()
print(finishing-starting)

## @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@  @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ ##
## Extract the number of common genes by subgroups and consider as true subgroups those with more common genes than the expected by chance ##
## @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@  @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ ##
## estimated duration: 50.48 mins
info<-read.csv2("Patients_information_table.txt",stringsAsFactors = F,sep="\t")
nocalc<-as.matrix(cbind(info[grep("Hutchin",info[,1]),c(1,2)],paste(gsub("_.+","",info[grep("Hutchin",info[,1]),1]),".",info[grep("Hutchin",info[,1]),2],sep="")))
starting<-Sys.time()
## generate a list of the patients and the subgroup they belong to
file_clusters<-list.files("Data/Genes_based_clustering/Clusters/")  ;  patients_subgroups<-c()
for(a in 1:length(file_clusters)){
  load(paste("Data/Genes_based_clustering/Clusters/",file_clusters[a],sep=""))
  patients_subgroups<-rbind(patients_subgroups,cbind(names(clusters$cluster),as.character(clusters$cluster),
                                                     paste(gsub("_.+","",names(clusters$cluster)),as.character(clusters$cluster),sep=".")))
}
colnames(nocalc)<-c("Patient","Cluster_number","Cluster_info")
colnames(patients_subgroups)<-c("Patient","Cluster_number","Cluster_info")
patients_subgroups<-rbind(patients_subgroups,nocalc)
## indicate the subgroups that are composed by only one patient
size1<-names(which(table(patients_subgroups[,3])==1))
for(a in 1:length(size1)){
  patients_subgroups[which(patients_subgroups[,3]==size1[a]),3]<-paste(strsplit(size1[a],".",fixed=T)[[1]][1],
                                                                       "^^Cluster_of_size_1_",strsplit(size1[a],".",fixed=T)[[1]][2],sep="")
}
colnames(patients_subgroups)<-c("Patient","Cluster_number","Cluster_info")
save(patients_subgroups,file="Data/First_subgroups.Rdata")
subgrupos<-unique(patients_subgroups[,3])
subgrupos<-subgrupos[-grep("^^",subgrupos,fixed=T)]
todos<-patients_subgroups[,c(1,3)]
todos<-todos[-grep("^^",todos[,2],fixed=T),]
## load the top 500 sDEGs in each patient
load("Data/Top_500_sDEGs_names.Rdata")
print("We start shuffling")
totales<-c()
listorras<-list()
for(z in 1:1000){
  genecom<-c()
  genescomunes<-list()
  for(a in 1:length(subgrupos)){
    coge<-subgrupos[a]
    ## randomly assign patients from the same disease to subgroups
    tos<-todos[grep(strsplit(coge,".",fixed=T)[[1]][1],todos[,2]),]  ;  todos2<-cbind(tos[,1],sample(tos[,2]))
    pacientes<-todos2[which(todos2[,2]==coge),1]  ;  upl<-list();downl<-list()
    for(b in 1:length(pacientes)){upl[[pacientes[b]]]<-degs[pacientes[b]][[1]]$up  ;  downl[[pacientes[b]]]<-degs[pacientes[b]][[1]]$down}
    ups<-Reduce(intersect,upl)  ;  downs<-Reduce(intersect,downl)  ;  genescomunes[[subgrupos[a]]]$up<-ups  ;  genescomunes[[subgrupos[a]]]$down<-downs
    if(length(ups)>0 || length(downs)>0){genecom<-c(genecom,1)}
  }
  totales<-c(totales,length(genecom))
  print(paste(z,length(genecom),sep="  "))
  listorras[[paste("Random_",z,sep="")]]<-genescomunes
}
print("Finished!!")
save(listorras,file="Data/Subgroups_common_genes_shuffled_intra_disease.Rdata")
finishing<-Sys.time()
finishing-starting

## @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ ##
## Extract those diseases whose subgroups present less common genes than the expected by chance ##
## @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ ##
load("Data/Top_500_sDEGs_names.Rdata") # degs
load("Data/First_subgroups.Rdata") # patients_subgroups
info<-read.csv2("Patients_information_table.txt",stringsAsFactors = F,sep="\t",header=T)
culur<-info$Disease_group_color  ;  names(culur)<-info$Disease
sizes<-as.numeric(table(patients_subgroups[,3])[colnames(plotbox)])  ;  names(sizes)<-colnames(plotbox)
## extract the genes shared by each patient-subgroup
groups<-unique(patients_subgroups[,3])
groups<-groups[-grep("^^Cluster",groups,fixed=T)]
commongenes<-list()
for(a in 1:length(groups)){
  ichoose<-patients_subgroups[which(patients_subgroups[,3]==groups[a]),1]
  ups<-list();downs<-list()  ;  for(b in 1:length(ichoose)){ups[[ichoose[b]]]<-degs[[ichoose[b]]]$up  ;  downs[[ichoose[b]]]<-degs[[ichoose[b]]]$down}
  commongenes[[groups[a]]]$up<-Reduce(intersect,ups)  ;  commongenes[[groups[a]]]$down<-Reduce(intersect,downs)
}
save(commongenes,file="Data/Subgroups_genes.Rdata")

## @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ ##
## Extract those subgroups presenting less common genes than the expected by chance ##
## @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ ##
## compare the genes shared by real subgroups vs. the shuffled ones
load("Data/Subgroups_genes.Rdata") # commongenes
load("Data/Subgroups_common_genes_shuffled_intra_disease.Rdata") # listorras
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

## @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @ @ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ ##
## We save the results with patients and their corresponding subgroup to generate the Patient_information... table ##
## @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @ @ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ ##
## select those clusters with less shared genes than the expected by chance and add a NS label
nosub<-colnames(plotbox)[which(color!="red")]
for(a in 1:length(nosub)){
  patients_subgroups[grep(nosub[a],patients_subgroups[,3]),2]<-paste(patients_subgroups[grep(nosub[a],patients_subgroups[,3]),2],"NS",sep="")
  patients_subgroups[grep(nosub[a],patients_subgroups[,3]),3]<-paste(gsub("\\..+","",patients_subgroups[grep(nosub[a],patients_subgroups[,3]),3]),".",
                                                                     patients_subgroups[grep(nosub[a],patients_subgroups[,3]),2],sep="")
}
clus<-patients_subgroups[,2]
names(clus)<-patients_subgroups[,1]
clusinf<-patients_subgroups[,3]
names(clusinf)<-patients_subgroups[,1]
## introduce the new subgroups into the patient_information table
info<-read.csv2("Patients_information_table.txt",stringsAsFactors = F,sep="\t")
if(length(grep("cluster_specific",colnames(info)))==0){
  info<-cbind(info,as.character(clus[info[,1]]),as.character(clusinf[info[,1]]))
  colnames(info)[c(14,15)]<-c("Cluster_number_cluster_specific","Cluster_info_cluster_specific")
  info<-info[,c(1,14,3:7,15,9:13,2,8)]
  write.table(info,"Patients_information_table.txt",row.names=F,quote=F,sep="\t")
}

