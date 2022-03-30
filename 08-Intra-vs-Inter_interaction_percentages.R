load("Data/Patient_patient_similarity_network.Rdata")

nit<-fisher_table[which(fisher_table[,3]!=0),]
info<-read.csv2("Patients_information_table.txt",stringsAsFactors = F,sep="\t",header=T)

pats<-info$Patient  ;  subs<-info$Cluster_info_cluster_specific  ;  names(subs)<-info$Patient
enfermedades<-c()  ;  subgrupos<-c()
for(a in 1:length(pats)){
  if(length(which(info$Disease==gsub("_.+","",pats[a])))>1){
    op1<-c();op2<-c()
    if(length(which(nit[,1]==pats[a]))>0){op1<-nit[which(nit[,1]==pats[a]),c(2,3)]  ;  colnames(op1)<-c("Pat","Inter")}
    if(length(which(nit[,2]==pats[a]))>0){op2<-nit[which(nit[,2]==pats[a]),c(1,3)]  ;  colnames(op2)<-c("Pat","Inter")}
    op<-rbind(op1,op2)
    intra<-as.character(op[grep(paste(gsub("_.+","",pats[a]),"_",sep=""),op[,1]),2])  ;  intraper<-length(which(intra==1))/(length(which(info$Disease==gsub("_.+","",pats[a])))-1)
    inter<-as.character(op[-grep(paste(gsub("_.+","",pats[a]),"_",sep=""),op[,1]),2])  ;  interper<-length(which(inter==1))/length(which(info$Disease!=gsub("_.+","",pats[a])))
    enfermedades<-t(c(pats[a],intraper,interper))  ;  write.table(enfermedades,"Data/ROC_diseases.txt",append=T,quote=F,sep="\t",row.names=F,col.names = F)
    if(length(which(info$Cluster_info_cluster_specific==as.character(subs[pats[a]])))>1){
      ops<-cbind(as.character(subs[op[,1]]),as.character(op[,2]))
      intra<-as.character(ops[which(ops[,1]==as.character(subs[pats[a]])),2])
      intraper2<-length(which(intra==1))/(length(which(info$Cluster_info_cluster_specific==as.character(subs[pats[a]])))-1)
      inter<-as.character(ops[which(ops[,1]!=as.character(subs[pats[a]])),2])
      interper2<-length(which(inter==1))/length(which(info$Cluster_info_cluster_specific!=as.character(subs[pats[a]])))
      subgrupos<-t(c(pats[a],intraper2,interper2))  ;  write.table(subgrupos,"Data/ROC_subgroups.txt",append=T,quote=F,sep="\t",row.names=F,col.names = F)
    }
  }
}

## @@ @@ @@ @@ @@ ##
## Generate plots ##
## @@ @@ @@ @@ @@ ##

val<-read.csv2("Data/ROC_subgroups.txt",stringsAsFactors = F,sep="\t",header=F)  ;  valores2<-val[,2:3]  ;  rownames(valores2)<-val[,1]
val<-read.csv2("Data/ROC_diseases.txt",stringsAsFactors = F,sep="\t",header=F)  ;  valores1<-val[,2:3]  ;  rownames(valores1)<-val[,1]

## Diseases ##
## @@ @@ @@ ##
valores_juntos<-c()
for(a in 1:length(valores1[,1])){
  if(valores1[a,1]<=1 && valores1[a,1]>0.9){valores_juntos<-rbind(valores_juntos,c(1,valores1[a,2]))}
  if(valores1[a,1]<=0.9 && valores1[a,1]>0.8){valores_juntos<-rbind(valores_juntos,c(0.9,valores1[a,2]))}
  if(valores1[a,1]<=0.8 && valores1[a,1]>0.7){valores_juntos<-rbind(valores_juntos,c(0.8,valores1[a,2]))}
  if(valores1[a,1]<=0.7 && valores1[a,1]>0.6){valores_juntos<-rbind(valores_juntos,c(0.7,valores1[a,2]))}
  if(valores1[a,1]<=0.6 && valores1[a,1]>0.5){valores_juntos<-rbind(valores_juntos,c(0.6,valores1[a,2]))}
  if(valores1[a,1]<=0.5 && valores1[a,1]>0.4){valores_juntos<-rbind(valores_juntos,c(0.5,valores1[a,2]))}
  if(valores1[a,1]<=0.4 && valores1[a,1]>0.3){valores_juntos<-rbind(valores_juntos,c(0.4,valores1[a,2]))}
  if(valores1[a,1]<=0.3 && valores1[a,1]>0.2){valores_juntos<-rbind(valores_juntos,c(0.3,valores1[a,2]))}
  if(valores1[a,1]<=0.2 && valores1[a,1]>0.1){valores_juntos<-rbind(valores_juntos,c(0.2,valores1[a,2]))}
  if(valores1[a,1]<=0.1 && valores1[a,1]>0){valores_juntos<-rbind(valores_juntos,c(0.1,valores1[a,2]))}
}
Xs<-c()
for(a in 1:length(valores_juntos[,1])){
  if(valores_juntos[a,1]==1){Xs<-c(Xs,10)}
  if(valores_juntos[a,1]==0.9){Xs<-c(Xs,9)}
  if(valores_juntos[a,1]==0.8){Xs<-c(Xs,8)}
  if(valores_juntos[a,1]==0.7){Xs<-c(Xs,7)}
  if(valores_juntos[a,1]==0.6){Xs<-c(Xs,6)}
  if(valores_juntos[a,1]==0.5){Xs<-c(Xs,5)}
  if(valores_juntos[a,1]==0.4){Xs<-c(Xs,4)}
  if(valores_juntos[a,1]==0.3){Xs<-c(Xs,3)}
  if(valores_juntos[a,1]==0.2){Xs<-c(Xs,2)}
  if(valores_juntos[a,1]==0.1){Xs<-c(Xs,1)}
}
limites <- cbind(seq(0,1,by = 0.1)[1:10], seq(0,1,by = 0.1)[2:11])  ;  Ys <- as.integer(NA)
for (i in 1:nrow(limites)) {m <- cbind(valores_juntos[,2] > limites[i, 1], valores_juntos[,2] <= limites[i, 2])  ;  Ys[apply(m, 1, all)] <- i}  ;  Ys[valores_juntos[,2]  == 0] <- 1
tabla <- table(Xs, Ys)

if(max(as.numeric(rownames(tabla)))<10){
  if(length(tabla[,1])<10){tabla<-rbind(tabla,matrix(ncol=length(tabla[1,]),nrow=10-length(tabla[,1]),0))}
  if(length(tabla[1,])<10){tabla<-cbind(tabla,matrix(ncol=10-length(tabla[1,]),nrow=length(tabla[,1]),0))}
}
if(max(as.numeric(rownames(tabla)))==10){
  if(length(tabla[,1])<10){tabla<-rbind(matrix(ncol=length(tabla[1,]),nrow=10-length(tabla[,1]),0),tabla)}
  if(length(tabla[1,])<10){tabla<-cbind(tabla,matrix(ncol=10-length(tabla[1,]),nrow=length(tabla[,1]),0))}
}
tabla2<-tabla[10:1,]

## Plot the obtained results ##
palette.breaks=c(seq(min(tabla2),-1,length=100),seq(-0.9999999,0.9999999,length=100),seq(1,max(tabla2),length=100))
color.palette<-colorRampPalette(c("#FF0000","#FFFFFF","#7d7e8b"))(length(palette.breaks)-1)
dev.off()
heatmap.2(tabla2,scale="none",dendrogram = "none",trace="none",col=color.palette,breaks=palette.breaks,
          Colv=FALSE,Rowv=FALSE,cexRow = 1,cexCol = 1,"key"=FALSE,margins = c(10,15),
          lmat=rbind(c(0,3,0),c(2,1,0),c(0,4,0)),lhei=c(0.5,8,0),lwid=c(0.5,8,0),cellnote=tabla2,notecex = 1.5,notecol="black")
pdf(file="Plots/ROCs/Diseases.pdf")
heatmap.2(tabla2,scale="none",dendrogram = "none",trace="none",col=color.palette,breaks=palette.breaks,
          Colv=FALSE,Rowv=FALSE,cexRow = 1,cexCol = 1,"key"=FALSE,margins = c(10,15),
          lmat=rbind(c(0,3,0),c(2,1,0),c(0,4,0)),lhei=c(0.5,8,0),lwid=c(0.5,8,0),cellnote=tabla2,notecex = 1.5,notecol="black")
dev.off()

## Clusters ##
## @@ @@ @@ ##
valores_juntos<-c()
for(a in 1:length(valores2[,1])){
  if(valores2[a,1]<=1 && valores2[a,1]>0.9){valores_juntos<-rbind(valores_juntos,c(1,valores2[a,2]))}
  if(valores2[a,1]<=0.9 && valores2[a,1]>0.8){valores_juntos<-rbind(valores_juntos,c(0.9,valores2[a,2]))}
  if(valores2[a,1]<=0.8 && valores2[a,1]>0.7){valores_juntos<-rbind(valores_juntos,c(0.8,valores2[a,2]))}
  if(valores2[a,1]<=0.7 && valores2[a,1]>0.6){valores_juntos<-rbind(valores_juntos,c(0.7,valores2[a,2]))}
  if(valores2[a,1]<=0.6 && valores2[a,1]>0.5){valores_juntos<-rbind(valores_juntos,c(0.6,valores2[a,2]))}
  if(valores2[a,1]<=0.5 && valores2[a,1]>0.4){valores_juntos<-rbind(valores_juntos,c(0.5,valores2[a,2]))}
  if(valores2[a,1]<=0.4 && valores2[a,1]>0.3){valores_juntos<-rbind(valores_juntos,c(0.4,valores2[a,2]))}
  if(valores2[a,1]<=0.3 && valores2[a,1]>0.2){valores_juntos<-rbind(valores_juntos,c(0.3,valores2[a,2]))}
  if(valores2[a,1]<=0.2 && valores2[a,1]>0.1){valores_juntos<-rbind(valores_juntos,c(0.2,valores2[a,2]))}
  if(valores2[a,1]<=0.1 && valores2[a,1]>0){valores_juntos<-rbind(valores_juntos,c(0.1,valores2[a,2]))}
}
Xs<-c()
for(a in 1:length(valores_juntos[,1])){
  if(valores_juntos[a,1]==1){Xs<-c(Xs,10)}
  if(valores_juntos[a,1]==0.9){Xs<-c(Xs,9)}
  if(valores_juntos[a,1]==0.8){Xs<-c(Xs,8)}
  if(valores_juntos[a,1]==0.7){Xs<-c(Xs,7)}
  if(valores_juntos[a,1]==0.6){Xs<-c(Xs,6)}
  if(valores_juntos[a,1]==0.5){Xs<-c(Xs,5)}
  if(valores_juntos[a,1]==0.4){Xs<-c(Xs,4)}
  if(valores_juntos[a,1]==0.3){Xs<-c(Xs,3)}
  if(valores_juntos[a,1]==0.2){Xs<-c(Xs,2)}
  if(valores_juntos[a,1]==0.1){Xs<-c(Xs,1)}
}
limites <- cbind(seq(0,1,by = 0.1)[1:10], seq(0,1,by = 0.1)[2:11])  ;  Ys <- as.integer(NA)
for (i in 1:nrow(limites)) {m <- cbind(valores_juntos[,2] > limites[i, 1], valores_juntos[,2] <= limites[i, 2])  ;  Ys[apply(m, 1, all)] <- i}  ;  Ys[valores_juntos[,2]  == 0] <- 1
tabla <- table(Xs, Ys)
if(max(as.numeric(rownames(tabla)))<10){
  if(length(tabla[,1])<10){tabla<-rbind(tabla,matrix(ncol=length(tabla[1,]),nrow=10-length(tabla[,1]),0))}
  if(length(tabla[1,])<10){tabla<-cbind(tabla,matrix(ncol=10-length(tabla[1,]),nrow=length(tabla[,1]),0))}
}
if(max(as.numeric(rownames(tabla)))==10){
  if(length(tabla[,1])<10){tabla<-rbind(matrix(ncol=length(tabla[1,]),nrow=10-length(tabla[,1]),0),tabla)}
  if(length(tabla[1,])<10){tabla<-cbind(tabla,matrix(ncol=10-length(tabla[1,]),nrow=length(tabla[,1]),0))}
}
tabla2<-tabla[10:1,]
## Plot the obtained results ##
palette.breaks=c(seq(min(tabla2),-1,length=100),seq(-0.9999999,0.9999999,length=100),seq(1,max(tabla2),length=100))
color.palette<-colorRampPalette(c("#1022ea","#FFFFFF","#7d7e8b"))(length(palette.breaks)-1)
dev.off()
heatmap.2(tabla2,scale="none",dendrogram = "none",trace="none",col=color.palette,breaks=palette.breaks,
          Colv=FALSE,Rowv=FALSE,cexRow = 1,cexCol = 1,"key"=FALSE,margins = c(10,15),
          lmat=rbind(c(0,3,0),c(2,1,0),c(0,4,0)),lhei=c(0.5,8,0),lwid=c(0.5,8,0),cellnote=tabla2,notecex = 1.5,notecol="black")
pdf(file="Plots/ROCs/Subgroups.pdf")
heatmap.2(tabla2,scale="none",dendrogram = "none",trace="none",col=color.palette,breaks=palette.breaks,
          Colv=FALSE,Rowv=FALSE,cexRow = 1,cexCol = 1,"key"=FALSE,margins = c(10,15),
          lmat=rbind(c(0,3,0),c(2,1,0),c(0,4,0)),lhei=c(0.5,8,0),lwid=c(0.5,8,0),cellnote=tabla2,notecex = 1.5,notecol="black")
dev.off()



## Diseases Mental ##
## @ @@ @@ @@ @@ @ ##
grupo<-info$Disease_group_name  ;  names(grupo)<-info$Patient
val<-read.csv2("Data/ROC_diseases.txt",stringsAsFactors = F,sep="\t",header=F)  ;  valores1<-val[,2:3]
cual1<-as.character(grupo[val[,1]])
valores1<-valores1[c(grep("Mental",cual1),grep("nervous",cual1)),]
valores_juntos<-c()
for(a in 1:length(valores1[,1])){
  if(valores1[a,1]<=1 && valores1[a,1]>0.9){valores_juntos<-rbind(valores_juntos,c(1,valores1[a,2]))}
  if(valores1[a,1]<=0.9 && valores1[a,1]>0.8){valores_juntos<-rbind(valores_juntos,c(0.9,valores1[a,2]))}
  if(valores1[a,1]<=0.8 && valores1[a,1]>0.7){valores_juntos<-rbind(valores_juntos,c(0.8,valores1[a,2]))}
  if(valores1[a,1]<=0.7 && valores1[a,1]>0.6){valores_juntos<-rbind(valores_juntos,c(0.7,valores1[a,2]))}
  if(valores1[a,1]<=0.6 && valores1[a,1]>0.5){valores_juntos<-rbind(valores_juntos,c(0.6,valores1[a,2]))}
  if(valores1[a,1]<=0.5 && valores1[a,1]>0.4){valores_juntos<-rbind(valores_juntos,c(0.5,valores1[a,2]))}
  if(valores1[a,1]<=0.4 && valores1[a,1]>0.3){valores_juntos<-rbind(valores_juntos,c(0.4,valores1[a,2]))}
  if(valores1[a,1]<=0.3 && valores1[a,1]>0.2){valores_juntos<-rbind(valores_juntos,c(0.3,valores1[a,2]))}
  if(valores1[a,1]<=0.2 && valores1[a,1]>0.1){valores_juntos<-rbind(valores_juntos,c(0.2,valores1[a,2]))}
  if(valores1[a,1]<=0.1 && valores1[a,1]>0){valores_juntos<-rbind(valores_juntos,c(0.1,valores1[a,2]))}
}
Xs<-c()
for(a in 1:length(valores_juntos[,1])){
  if(valores_juntos[a,1]==1){Xs<-c(Xs,10)}
  if(valores_juntos[a,1]==0.9){Xs<-c(Xs,9)}
  if(valores_juntos[a,1]==0.8){Xs<-c(Xs,8)}
  if(valores_juntos[a,1]==0.7){Xs<-c(Xs,7)}
  if(valores_juntos[a,1]==0.6){Xs<-c(Xs,6)}
  if(valores_juntos[a,1]==0.5){Xs<-c(Xs,5)}
  if(valores_juntos[a,1]==0.4){Xs<-c(Xs,4)}
  if(valores_juntos[a,1]==0.3){Xs<-c(Xs,3)}
  if(valores_juntos[a,1]==0.2){Xs<-c(Xs,2)}
  if(valores_juntos[a,1]==0.1){Xs<-c(Xs,1)}
}
limites <- cbind(seq(0,1,by = 0.1)[1:10], seq(0,1,by = 0.1)[2:11])  ;  Ys <- as.integer(NA)
for (i in 1:nrow(limites)) {m <- cbind(valores_juntos[,2] > limites[i, 1], valores_juntos[,2] <= limites[i, 2])  ;  Ys[apply(m, 1, all)] <- i}  ;  Ys[valores_juntos[,2]  == 0] <- 1
tabla <- table(Xs, Ys)

if(max(as.numeric(rownames(tabla)))<10){
  if(length(tabla[,1])<10){tabla<-rbind(tabla,matrix(ncol=length(tabla[1,]),nrow=10-length(tabla[,1]),0))}
  if(length(tabla[1,])<10){tabla<-cbind(tabla,matrix(ncol=10-length(tabla[1,]),nrow=length(tabla[,1]),0))}
}
if(max(as.numeric(rownames(tabla)))==10){
  if(length(tabla[,1])<10){tabla<-rbind(matrix(ncol=length(tabla[1,]),nrow=10-length(tabla[,1]),0),tabla)}
  if(length(tabla[1,])<10){tabla<-cbind(tabla,matrix(ncol=10-length(tabla[1,]),nrow=length(tabla[,1]),0))}
}

tabla2<-tabla[10:1,]
## Plot the obtained results ##
palette.breaks=c(seq(min(tabla2),-1,length=100),seq(-0.9999999,0.9999999,length=100),seq(1,max(tabla2),length=100))
color.palette<-colorRampPalette(c("#FF0000","#FFFFFF","#7d7e8b"))(length(palette.breaks)-1)
dev.off()
heatmap.2(tabla2,scale="none",dendrogram = "none",trace="none",col=color.palette,breaks=palette.breaks,
          Colv=FALSE,Rowv=FALSE,cexRow = 1,cexCol = 1,"key"=FALSE,margins = c(10,15),
          lmat=rbind(c(0,3,0),c(2,1,0),c(0,4,0)),lhei=c(0.5,8,0),lwid=c(0.5,8,0),cellnote=tabla2,notecex = 1.5,notecol="black")
pdf(file="Plots/ROCs/Diseases_mental.pdf")
heatmap.2(tabla2,scale="none",dendrogram = "none",trace="none",col=color.palette,breaks=palette.breaks,
          Colv=FALSE,Rowv=FALSE,cexRow = 1,cexCol = 1,"key"=FALSE,margins = c(10,15),
          lmat=rbind(c(0,3,0),c(2,1,0),c(0,4,0)),lhei=c(0.5,8,0),lwid=c(0.5,8,0),cellnote=tabla2,notecex = 1.5,notecol="black")
dev.off()


## Diseases Neoplasm ##
## @@ @@ @@ @@ @@ @@ ##
grupo<-info$Disease_group_name  ;  names(grupo)<-info$Patient
val<-read.csv2("Data/ROC_diseases.txt",stringsAsFactors = F,sep="\t",header=F)  ;  valores1<-val[,2:3]
cual1<-as.character(grupo[val[,1]])
valores1<-valores1[grep("Neoplas",cual1),]
valores_juntos<-c()
for(a in 1:length(valores1[,1])){
  if(valores1[a,1]<=1 && valores1[a,1]>0.9){valores_juntos<-rbind(valores_juntos,c(1,valores1[a,2]))}
  if(valores1[a,1]<=0.9 && valores1[a,1]>0.8){valores_juntos<-rbind(valores_juntos,c(0.9,valores1[a,2]))}
  if(valores1[a,1]<=0.8 && valores1[a,1]>0.7){valores_juntos<-rbind(valores_juntos,c(0.8,valores1[a,2]))}
  if(valores1[a,1]<=0.7 && valores1[a,1]>0.6){valores_juntos<-rbind(valores_juntos,c(0.7,valores1[a,2]))}
  if(valores1[a,1]<=0.6 && valores1[a,1]>0.5){valores_juntos<-rbind(valores_juntos,c(0.6,valores1[a,2]))}
  if(valores1[a,1]<=0.5 && valores1[a,1]>0.4){valores_juntos<-rbind(valores_juntos,c(0.5,valores1[a,2]))}
  if(valores1[a,1]<=0.4 && valores1[a,1]>0.3){valores_juntos<-rbind(valores_juntos,c(0.4,valores1[a,2]))}
  if(valores1[a,1]<=0.3 && valores1[a,1]>0.2){valores_juntos<-rbind(valores_juntos,c(0.3,valores1[a,2]))}
  if(valores1[a,1]<=0.2 && valores1[a,1]>0.1){valores_juntos<-rbind(valores_juntos,c(0.2,valores1[a,2]))}
  if(valores1[a,1]<=0.1 && valores1[a,1]>0){valores_juntos<-rbind(valores_juntos,c(0.1,valores1[a,2]))}
}
Xs<-c()
for(a in 1:length(valores_juntos[,1])){
  if(valores_juntos[a,1]==1){Xs<-c(Xs,10)}
  if(valores_juntos[a,1]==0.9){Xs<-c(Xs,9)}
  if(valores_juntos[a,1]==0.8){Xs<-c(Xs,8)}
  if(valores_juntos[a,1]==0.7){Xs<-c(Xs,7)}
  if(valores_juntos[a,1]==0.6){Xs<-c(Xs,6)}
  if(valores_juntos[a,1]==0.5){Xs<-c(Xs,5)}
  if(valores_juntos[a,1]==0.4){Xs<-c(Xs,4)}
  if(valores_juntos[a,1]==0.3){Xs<-c(Xs,3)}
  if(valores_juntos[a,1]==0.2){Xs<-c(Xs,2)}
  if(valores_juntos[a,1]==0.1){Xs<-c(Xs,1)}
}
limites <- cbind(seq(0,1,by = 0.1)[1:10], seq(0,1,by = 0.1)[2:11])  ;  Ys <- as.integer(NA)
for (i in 1:nrow(limites)) {m <- cbind(valores_juntos[,2] > limites[i, 1], valores_juntos[,2] <= limites[i, 2])  ;  Ys[apply(m, 1, all)] <- i}  ;  Ys[valores_juntos[,2]  == 0] <- 1
tabla <- table(Xs, Ys)
if(max(as.numeric(rownames(tabla)))<10){
  if(length(tabla[,1])<10){tabla<-rbind(tabla,matrix(ncol=length(tabla[1,]),nrow=10-length(tabla[,1]),0))}
  if(length(tabla[1,])<10){tabla<-cbind(tabla,matrix(ncol=10-length(tabla[1,]),nrow=length(tabla[,1]),0))}
}
if(max(as.numeric(rownames(tabla)))==10){
  if(length(tabla[,1])<10){tabla<-rbind(matrix(ncol=length(tabla[1,]),nrow=10-length(tabla[,1]),0),tabla)}
  if(length(tabla[1,])<10){tabla<-cbind(tabla,matrix(ncol=10-length(tabla[1,]),nrow=length(tabla[,1]),0))}
}
tabla2<-tabla[10:1,]
## Plot the obtained results ##
palette.breaks=c(seq(min(tabla2),-1,length=100),seq(-0.9999999,0.9999999,length=100),seq(1,max(tabla2),length=100))
color.palette<-colorRampPalette(c("#FF0000","#FFFFFF","#7d7e8b"))(length(palette.breaks)-1)
dev.off()
heatmap.2(tabla2,scale="none",dendrogram = "none",trace="none",col=color.palette,breaks=palette.breaks,
          Colv=FALSE,Rowv=FALSE,cexRow = 1,cexCol = 1,"key"=FALSE,margins = c(10,15),
          lmat=rbind(c(0,3,0),c(2,1,0),c(0,4,0)),lhei=c(0.5,8,0),lwid=c(0.5,8,0),cellnote=tabla2,notecex = 1.5,notecol="black")
pdf(file="Plots/ROCs/Diseases_neoplasms.pdf")
heatmap.2(tabla2,scale="none",dendrogram = "none",trace="none",col=color.palette,breaks=palette.breaks,
          Colv=FALSE,Rowv=FALSE,cexRow = 1,cexCol = 1,"key"=FALSE,margins = c(10,15),
          lmat=rbind(c(0,3,0),c(2,1,0),c(0,4,0)),lhei=c(0.5,8,0),lwid=c(0.5,8,0),cellnote=tabla2,notecex = 1.5,notecol="black")
dev.off()






## Clusters mental ##
## @@ @@ @ @ @@ @@ ##
grupo<-info$Disease_group_name  ;  names(grupo)<-info$Patient
val<-read.csv2("Data/ROC_subgroups.txt",stringsAsFactors = F,sep="\t",header=F)  ;  valores2<-val[,2:3]
cual2<-as.character(grupo[val[,1]])
valores2<-valores2[c(grep("Mental",cual2),grep("nervous",cual2)),]
valores_juntos<-c()
for(a in 1:length(valores2[,1])){
  if(valores2[a,1]<=1 && valores2[a,1]>0.9){valores_juntos<-rbind(valores_juntos,c(1,valores2[a,2]))}
  if(valores2[a,1]<=0.9 && valores2[a,1]>0.8){valores_juntos<-rbind(valores_juntos,c(0.9,valores2[a,2]))}
  if(valores2[a,1]<=0.8 && valores2[a,1]>0.7){valores_juntos<-rbind(valores_juntos,c(0.8,valores2[a,2]))}
  if(valores2[a,1]<=0.7 && valores2[a,1]>0.6){valores_juntos<-rbind(valores_juntos,c(0.7,valores2[a,2]))}
  if(valores2[a,1]<=0.6 && valores2[a,1]>0.5){valores_juntos<-rbind(valores_juntos,c(0.6,valores2[a,2]))}
  if(valores2[a,1]<=0.5 && valores2[a,1]>0.4){valores_juntos<-rbind(valores_juntos,c(0.5,valores2[a,2]))}
  if(valores2[a,1]<=0.4 && valores2[a,1]>0.3){valores_juntos<-rbind(valores_juntos,c(0.4,valores2[a,2]))}
  if(valores2[a,1]<=0.3 && valores2[a,1]>0.2){valores_juntos<-rbind(valores_juntos,c(0.3,valores2[a,2]))}
  if(valores2[a,1]<=0.2 && valores2[a,1]>0.1){valores_juntos<-rbind(valores_juntos,c(0.2,valores2[a,2]))}
  if(valores2[a,1]<=0.1 && valores2[a,1]>0){valores_juntos<-rbind(valores_juntos,c(0.1,valores2[a,2]))}
}
Xs<-c()
for(a in 1:length(valores_juntos[,1])){
  if(valores_juntos[a,1]==1){Xs<-c(Xs,10)}
  if(valores_juntos[a,1]==0.9){Xs<-c(Xs,9)}
  if(valores_juntos[a,1]==0.8){Xs<-c(Xs,8)}
  if(valores_juntos[a,1]==0.7){Xs<-c(Xs,7)}
  if(valores_juntos[a,1]==0.6){Xs<-c(Xs,6)}
  if(valores_juntos[a,1]==0.5){Xs<-c(Xs,5)}
  if(valores_juntos[a,1]==0.4){Xs<-c(Xs,4)}
  if(valores_juntos[a,1]==0.3){Xs<-c(Xs,3)}
  if(valores_juntos[a,1]==0.2){Xs<-c(Xs,2)}
  if(valores_juntos[a,1]==0.1){Xs<-c(Xs,1)}
}
limites <- cbind(seq(0,1,by = 0.1)[1:10], seq(0,1,by = 0.1)[2:11])  ;  Ys <- as.integer(NA)
for (i in 1:nrow(limites)) {m <- cbind(valores_juntos[,2] > limites[i, 1], valores_juntos[,2] <= limites[i, 2])  ;  Ys[apply(m, 1, all)] <- i}  ;  Ys[valores_juntos[,2]  == 0] <- 1
tabla <- table(Xs, Ys)
if(max(as.numeric(rownames(tabla)))<10){
  if(length(tabla[,1])<10){tabla<-rbind(tabla,matrix(ncol=length(tabla[1,]),nrow=10-length(tabla[,1]),0))}
  if(length(tabla[1,])<10){tabla<-cbind(tabla,matrix(ncol=10-length(tabla[1,]),nrow=length(tabla[,1]),0))}
}
if(max(as.numeric(rownames(tabla)))==10){
  if(length(tabla[,1])<10){tabla<-rbind(matrix(ncol=length(tabla[1,]),nrow=10-length(tabla[,1]),0),tabla)}
  if(length(tabla[1,])<10){tabla<-cbind(tabla,matrix(ncol=10-length(tabla[1,]),nrow=length(tabla[,1]),0))}
}
tabla2<-tabla[10:1,]
## Plot the obtained results ##
palette.breaks=c(seq(min(tabla2),-1,length=100),seq(-0.9999999,0.9999999,length=100),seq(1,max(tabla2),length=100))
color.palette<-colorRampPalette(c("#1022ea","#FFFFFF","#7d7e8b"))(length(palette.breaks)-1)
dev.off()
heatmap.2(tabla2,scale="none",dendrogram = "none",trace="none",col=color.palette,breaks=palette.breaks,
          Colv=FALSE,Rowv=FALSE,cexRow = 1,cexCol = 1,"key"=FALSE,margins = c(10,15),
          lmat=rbind(c(0,3,0),c(2,1,0),c(0,4,0)),lhei=c(0.5,8,0),lwid=c(0.5,8,0),cellnote=tabla2,notecex = 1.5,notecol="black")
pdf(file="Plots/ROCs/Subgroups_mental.pdf")
heatmap.2(tabla2,scale="none",dendrogram = "none",trace="none",col=color.palette,breaks=palette.breaks,
          Colv=FALSE,Rowv=FALSE,cexRow = 1,cexCol = 1,"key"=FALSE,margins = c(10,15),
          lmat=rbind(c(0,3,0),c(2,1,0),c(0,4,0)),lhei=c(0.5,8,0),lwid=c(0.5,8,0),cellnote=tabla2,notecex = 1.5,notecol="black")
dev.off()


## Clusters neoplasms ##
## @@ @@ @@  @@ @@ @@ ##
grupo<-info$Disease_group_name  ;  names(grupo)<-info$Patient
val<-read.csv2("Data/ROC_subgroups.txt",stringsAsFactors = F,sep="\t",header=F)  ;  valores2<-val[,2:3]
cual2<-as.character(grupo[val[,1]])
valores2<-valores2[grep("Neoplas",cual2),]
valores_juntos<-c()
for(a in 1:length(valores2[,1])){
  if(valores2[a,1]<=1 && valores2[a,1]>0.9){valores_juntos<-rbind(valores_juntos,c(1,valores2[a,2]))}
  if(valores2[a,1]<=0.9 && valores2[a,1]>0.8){valores_juntos<-rbind(valores_juntos,c(0.9,valores2[a,2]))}
  if(valores2[a,1]<=0.8 && valores2[a,1]>0.7){valores_juntos<-rbind(valores_juntos,c(0.8,valores2[a,2]))}
  if(valores2[a,1]<=0.7 && valores2[a,1]>0.6){valores_juntos<-rbind(valores_juntos,c(0.7,valores2[a,2]))}
  if(valores2[a,1]<=0.6 && valores2[a,1]>0.5){valores_juntos<-rbind(valores_juntos,c(0.6,valores2[a,2]))}
  if(valores2[a,1]<=0.5 && valores2[a,1]>0.4){valores_juntos<-rbind(valores_juntos,c(0.5,valores2[a,2]))}
  if(valores2[a,1]<=0.4 && valores2[a,1]>0.3){valores_juntos<-rbind(valores_juntos,c(0.4,valores2[a,2]))}
  if(valores2[a,1]<=0.3 && valores2[a,1]>0.2){valores_juntos<-rbind(valores_juntos,c(0.3,valores2[a,2]))}
  if(valores2[a,1]<=0.2 && valores2[a,1]>0.1){valores_juntos<-rbind(valores_juntos,c(0.2,valores2[a,2]))}
  if(valores2[a,1]<=0.1 && valores2[a,1]>0){valores_juntos<-rbind(valores_juntos,c(0.1,valores2[a,2]))}
}
Xs<-c()
for(a in 1:length(valores_juntos[,1])){
  if(valores_juntos[a,1]==1){Xs<-c(Xs,10)}
  if(valores_juntos[a,1]==0.9){Xs<-c(Xs,9)}
  if(valores_juntos[a,1]==0.8){Xs<-c(Xs,8)}
  if(valores_juntos[a,1]==0.7){Xs<-c(Xs,7)}
  if(valores_juntos[a,1]==0.6){Xs<-c(Xs,6)}
  if(valores_juntos[a,1]==0.5){Xs<-c(Xs,5)}
  if(valores_juntos[a,1]==0.4){Xs<-c(Xs,4)}
  if(valores_juntos[a,1]==0.3){Xs<-c(Xs,3)}
  if(valores_juntos[a,1]==0.2){Xs<-c(Xs,2)}
  if(valores_juntos[a,1]==0.1){Xs<-c(Xs,1)}
}
limites <- cbind(seq(0,1,by = 0.1)[1:10], seq(0,1,by = 0.1)[2:11])  ;  Ys <- as.integer(NA)
for (i in 1:nrow(limites)) {m <- cbind(valores_juntos[,2] > limites[i, 1], valores_juntos[,2] <= limites[i, 2])  ;  Ys[apply(m, 1, all)] <- i}  ;  Ys[valores_juntos[,2]  == 0] <- 1
tabla <- table(Xs, Ys)
if(max(as.numeric(rownames(tabla)))<10){
  if(length(tabla[,1])<10){tabla<-rbind(tabla,matrix(ncol=length(tabla[1,]),nrow=10-length(tabla[,1]),0))}
  if(length(tabla[1,])<10){tabla<-cbind(tabla,matrix(ncol=10-length(tabla[1,]),nrow=length(tabla[,1]),0))}
}
if(max(as.numeric(rownames(tabla)))==10){
  if(length(tabla[,1])<10){tabla<-rbind(matrix(ncol=length(tabla[1,]),nrow=10-length(tabla[,1]),0),tabla)}
  if(length(tabla[1,])<10){tabla<-cbind(tabla,matrix(ncol=10-length(tabla[1,]),nrow=length(tabla[,1]),0))}
}
tabla2<-tabla[10:1,]
## Plot the obtained results ##
palette.breaks=c(seq(min(tabla2),-1,length=100),seq(-0.9999999,0.9999999,length=100),seq(1,max(tabla2),length=100))
color.palette<-colorRampPalette(c("#1022ea","#FFFFFF","#7d7e8b"))(length(palette.breaks)-1)
dev.off()
heatmap.2(tabla2,scale="none",dendrogram = "none",trace="none",col=color.palette,breaks=palette.breaks,
          Colv=FALSE,Rowv=FALSE,cexRow = 1,cexCol = 1,"key"=FALSE,margins = c(10,15),
          lmat=rbind(c(0,3,0),c(2,1,0),c(0,4,0)),lhei=c(0.5,8,0),lwid=c(0.5,8,0),cellnote=tabla2,notecex = 1.5,notecol="black")
pdf(file="Plots/ROCs/Subgroups_neoplasms.pdf")
heatmap.2(tabla2,scale="none",dendrogram = "none",trace="none",col=color.palette,breaks=palette.breaks,
          Colv=FALSE,Rowv=FALSE,cexRow = 1,cexCol = 1,"key"=FALSE,margins = c(10,15),
          lmat=rbind(c(0,3,0),c(2,1,0),c(0,4,0)),lhei=c(0.5,8,0),lwid=c(0.5,8,0),cellnote=tabla2,notecex = 1.5,notecol="black")
dev.off()



