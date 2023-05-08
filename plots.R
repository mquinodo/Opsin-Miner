

#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

file1=args[1]
file2=args[2]
file3=args[3]
file4=args[4]
file5=args[5]
dir=args[6]

load(file5)

# reading files
data=read.table(file=file1,header=T,sep="\t",check.names = FALSE)
haplo=read.table(file=file2)
haploP=read.table(file=file3)
rareVariants=read.table(file=file4,header=T,row.names=NULL,sep="\t",colClasses = c("character"))

# processing variants
rareVariants=cbind(rareVariants,rareVariants[,1:3])
colnames(rareVariants)[9:11]=c("HGVS-g.","HGVS-c.","HGVS-p.")
for (i in 1:dim(rareVariants)[1]){
	sel=which(HGVS[,7]==rareVariants[i,3] & HGVS[,8]==rareVariants[i,4] & HGVS[,9]==rareVariants[i,5])
	if(length(sel)==1){
		rareVariants[i,9:11]=HGVS[sel,c(4,2,3)]
	} else {
		rareVariants[i,9:11]=c("NA","NA","Check c. and p. notations")
	}
}
rareVariants=rareVariants[which(rareVariants[,9]!="NA" | nchar(rareVariants[,4])>1 | nchar(rareVariants[,5])>1),]

# normalize by chrX and remultiply by mean
data2=data
for (i in 2:dim(data2)[2]){
	data2[3:26,i]=as.numeric(data[3:26,i])/as.numeric(data[2,i])
}
data3=data2[,2:dim(data2)[2]]
rownames(data3)=data2[,1]
for (i in 1:dim(data3)[2]){
	data3[,i]=as.numeric(data3[,i])*mean(as.numeric(data[2,2:dim(data)[2]]))
}

# extracing coverage info for various exons
LW=1:dim(data3)[2]
MWMW2=1:dim(data3)[2]
MW=1:dim(data3)[2]
MW2=1:dim(data3)[2]
exons<-matrix(nrow=12,ncol=dim(data3)[2])
rownames(exons)=c("LCR*","ALL\nexon 3","ALL\nexon 6","OPN1LW\nexon 1","OPN1LW\nexon 2","OPN1LW\nexon 4","OPN1LW\nexon 5","OPN1MW\nexon 5","OPN1MW2\nexon 5","OPN1MW/MW2\nexon 1","OPN1MW/MW2\nexon 2","OPN1MW/MW2\nexon 4")
namesexons=c("LCR","ALL-exon3","ALLexon-6","OPN1LW-exon1","OPN1LW-exon2","OPN1LW-exon4","OPN1LW-exon5","OPN1MW-exon5","OPN1MW2-exon5","OPN1MW/MW2-exon1","OPN1MW/MW2-exon2","OPN1MW/MW2-exon4")
exons2=exons

# mean per target without sample with less than 10% of average
me=data3[,1]
for (i in 1:dim(data3)[1]){
	me[i]=mean(as.numeric(data3[i,which(as.numeric(data3[i,])>0.1*mean(as.numeric(data3[i,])))]))
}

data4=data[,-1]
for (i in 1:dim(data3)[2]){
	LW[i]=sum(data3[9:16,i])/sum(me[9:16])
	MWMW2[i]=sum(data3[21:26,i])/sum(me[21:26])
	MW[i]=sum(data3[c(17,18),i])/sum(me[17:18])
	MW2[i]=sum(data3[c(19,20),i])/sum(me[19:20])

	exons[1,i]=max(as.numeric(data3[3,i])/me[3])
	exons[2,i]=max(as.numeric(data3[5,i])/me[5],as.numeric(data3[6,i])/me[6])
	exons[3,i]=max(as.numeric(data3[7,i])/me[7],as.numeric(data3[8,i])/me[8])
	exons[4,i]=max(as.numeric(data3[9,i])/me[9],as.numeric(data3[10,i])/me[10])
	exons[5,i]=max(as.numeric(data3[11,i])/me[11],as.numeric(data3[12,i])/me[12])
	exons[6,i]=max(as.numeric(data3[13,i])/me[13],as.numeric(data3[14,i])/me[14])
	exons[7,i]=max(as.numeric(data3[15,i])/me[15],as.numeric(data3[16,i])/me[16])
	exons[8,i]=max(as.numeric(data3[17,i])/me[17],as.numeric(data3[18,i])/me[18])
	exons[9,i]=max(as.numeric(data3[19,i])/me[19],as.numeric(data3[20,i])/me[20])
	exons[10,i]=max(as.numeric(data3[21,i])/me[21],as.numeric(data3[22,i])/me[22])
	exons[11,i]=max(as.numeric(data3[23,i])/me[23],as.numeric(data3[24,i])/me[24])
	exons[12,i]=max(as.numeric(data3[25,i])/me[25],as.numeric(data3[26,i])/me[26])

	exons2[1,i]=max(as.numeric(data4[3,i]))
	exons2[2,i]=max(as.numeric(data4[5,i]),as.numeric(data4[6,i]))
	exons2[3,i]=max(as.numeric(data4[7,i]),as.numeric(data4[8,i]))
	exons2[4,i]=max(as.numeric(data4[9,i]),as.numeric(data4[10,i]))
	exons2[5,i]=max(as.numeric(data4[11,i]),as.numeric(data4[12,i]))
	exons2[6,i]=max(as.numeric(data4[13,i]),as.numeric(data4[14,i]))
	exons2[7,i]=max(as.numeric(data4[15,i]),as.numeric(data4[16,i]))
	exons2[8,i]=max(as.numeric(data4[17,i]),as.numeric(data4[18,i]))
	exons2[9,i]=max(as.numeric(data4[19,i]),as.numeric(data4[20,i]))
	exons2[10,i]=max(as.numeric(data4[21,i]),as.numeric(data4[22,i]))
	exons2[11,i]=max(as.numeric(data4[23,i]),as.numeric(data4[24,i]))
	exons2[12,i]=max(as.numeric(data4[25,i]),as.numeric(data4[26,i]))
}

medi1=c(median(apply(data4[9:16,],2,sum)),median(apply(data4[21:26,],2,sum)),median(apply(data4[17:18,],2,sum)),median(apply(data4[19:20,],2,sum)))
medi2=c(median(as.numeric(data4[3,])),median(apply(data4[5:6,],2,max)),median(apply(data4[7:8,],2,max)),median(apply(data4[9:10,],2,max)),
	median(apply(data4[11:12,],2,max)),median(apply(data4[13:14,],2,max)),median(apply(data4[15:16,],2,max)),median(apply(data4[17:18,],2,max)),
	median(apply(data4[19:20,],2,max)),median(apply(data4[21:22,],2,max)),median(apply(data4[23:24,],2,max)),median(apply(data4[25:26,],2,max)))

exons=exons[c(1,4:12,2:3),]
exons2=exons2[c(1,4:12,2:3),]
medi2=medi2[c(1,4:12,2:3)]
namesexons=namesexons[c(1,4:12,2:3)]

# inferring sex
rat=as.numeric(data[2,2:dim(data)[2]])/as.numeric(data[1,2:dim(data)[2]])
a=seq(min(rat,na.rm=T),max(rat,na.rm=T),(max(rat,na.rm=T)-min(rat,na.rm=T))/100)
b<-matrix(10,ncol=101,nrow=101)
for (i in 1:length(a)){
	for (j in 1:length(a)){
		if(i<j){
			a1=rat[which(rat>(a[i]+a[j])/2)]
			a2=rat[which(rat<(a[i]+a[j])/2)]
			b[i,j]=sum(abs(a1-a[j])^2)+sum(abs(a2-a[i])^2)
		}
	}
}
c=which(b == min(b,na.rm=T), arr.ind = TRUE)
zm=rat
zf=rat
lim1=min(a[c[1]],a[c[2]])
lim2=max(a[c[1]],a[c[2]])
m=which(rat<(3*lim1+lim2)/4 & rat>(5*lim1-lim2)/4)
f=which(rat<(5*lim2-lim1)/4 & rat>(3*lim2+lim1)/4)
for (i in 1:length(rat)){
	zm[i]=(rat[i]-mean(rat[m]))/sd(rat[m])
	zf[i]=(rat[i]-mean(rat[f]))/sd(rat[f])
}
z1=zm+zf
sex=z1
sex[which(z1<(-3))]="Male"
sex[which(z1>(3))]="Female"
sex[which(z1<=(3) & z1>=(-3))]="Undetermined"

pdf(file=paste(dir,"/06_plots/","0.sex.pdf",sep=""),height=10,width=20)
plot(z1)
dev.off()

sexout=cbind(colnames(data3),sex,z1)
colnames(sexout)=c("ID","Inferred-sex","Z-score")
write.table(sexout,file=paste(dir,"/0.sex.tsv",sep=""),quote=F,row.names=F,sep="\t")

events<-matrix(ncol=9,nrow=0)
colnames(events)=c("#ID","Inferred-sex","Event","Pathogenicity","Details","Supporting-reads","WT-reads","Ratio","gnomAD-AF")

phenotype<-matrix(ncol=4,nrow=0)
colnames(phenotype)=c("#ID","Inferred-sex","Score","Phenotype")

# plots
for (i in 1:dim(data3)[2]){

	score1=0
	score2=0
	score3=0

	pdf(file=paste(dir,"/06_plots/",colnames(data3)[i],"-both.pdf",sep=""),height=10,width=20)
	par(mfrow=c(1,2),mgp=c(2.5,1,0),mar=c(8.1, 4.1, 4.1, 2.1))
	m=max(2,c(LW,MWMW2,MW,MW2)*1,na.rm=T)
	m=max(2.5,c(LW[i],MWMW2[i],MW[i],MW2[i])*1.2,na.rm=T)
	plot(1:4,ylim=c(0,m),xlim=c(0.5,6.8), xaxt = "n",xlab="",ylab="Coverage ratio vs. controls",cex=0,main=paste("OPN1LW, MW and MW2 analysis for ",colnames(data3)[i],sep=""),cex.main=2,cex.lab=1.4)
	par(mgp=c(2.5,3,0))

	for (j in 1:4){
		if(medi1[j]<20){
			rect(j-0.4,-0.05,j+0.4,m*1.02,border=NA,col="#ffcccc")
		}
		if(medi1[j]>=20 & medi1[j]<=100){
			rect(j-0.4,-0.05,j+0.4,m*1.02,border=NA,col="#fde7b5")
		}
	}

	for (j in 0:10){
		lines(c(0,4.4),c(j,j),lty=2,col="darkgrey")
	}
	lines(c(0,4.4),c(1,1))
	lines(c(0,4.4),c(0,0))
	cols=c(4,4,4,4)
	if(LW[i]<0.1){cols[1]=2}
	if(MWMW2[i]<0.1){cols[2]=2}
	if(MW[i]<0.1){cols[3]=2}
	if(MW2[i]<0.1){cols[4]=2}

	axis(1, at=c(1,2,3,4), labels=c(
		paste("OPN1LW\nexons 1,2,4,5\nratio = ",round(LW[i],digits=2),sep=""),
		paste("OPN1MW/MW2\nexons 1,2,4\nratio = ",round(MWMW2[i],digits=2),sep=""),
		paste("OPN1MW\nexon 5\nratio = ",round(MW[i],digits=2),sep=""),
		paste("OPN1MW2\nexon 5\nratio = ",round(MW2[i],digits=2),sep="")))
	n=length(LW)-1

	rect(0.7,0,1.3,LW[i],col=cols[1])
	rect(1.7,0,2.3,MWMW2[i],col=cols[2])
	rect(2.7,0,3.3,MW[i],col=cols[3])
	rect(3.7,0,4.3,MW2[i],col=cols[4])

	if(cols[1]!=4 & LW[i]<0.02){lines(c(0.7,1.3),c(LW[i],LW[i]),col=cols[1],lwd=6)}
	if(cols[2]!=4 & MWMW2[i]<0.02){lines(c(1.7,2.3),c(MWMW2[i],MWMW2[i]),col=cols[2],lwd=6)}
	if(cols[3]!=4 & MW[i]<0.02){lines(c(2.7,3.3),c(MW[i],MW[i]),col=cols[3],lwd=6)}
	if(cols[4]!=4 & MW2[i]<0.02){lines(c(3.7,4.3),c(MW2[i],MW2[i]),col=cols[4],lwd=6)}

	points(rep(1,n)+runif(n,-0.2,0.2),LW[-i],pch=16,cex=0.7)
	points(rep(2,n)+runif(n,-0.2,0.2),MWMW2[-i],pch=16,cex=0.8)
	points(rep(3,n)+runif(n,-0.2,0.2),MW[-i],pch=16,cex=0.8)
	points(rep(4,n)+runif(n,-0.2,0.2),MW2[-i],pch=16,cex=0.8)

	text(4.7,m*0.89,paste("Inferred sex: ",sex[i],sep=""),cex=1.2,adj=0)
	legend(4.7,m*0.85,c("No deletion","Deletion","Other individuals","Very low coverage (<10)","Low coverage (<20)"),col=c(4,2,1,"#ffcccc","#fde7b5"),pch=c(15,15,16,15,15),pt.cex=c(2,2,1,2,2),bg="white",y.intersp=1.5)

	if(LW[i]<0.1 & MWMW2[i]>=0.1){
		text(4.7,m*0.02,"Potential deletion\nof OPN1LW",adj=0,col=2,cex=1.2)
	}
	if(LW[i]>=0.1 & MWMW2[i]<0.1){
		text(4.7,m*0.02,"Potential deletion\nof OPN1MW/MW2",adj=0,col=2,cex=1.2)
	}
	if(LW[i]<0.1 & MWMW2[i]<0.1){
		text(4.7,m*0.05,"Potential deletion\nof OPN1LW and\nOPN1MW/MW2",adj=0,col=2,cex=1.2)
	}
	if(LW[i]>=0.1 & MWMW2[i]>=0.1 & MW[i]>=0.1 & MW2[i]>=0.1){
		text(4.7,m*0.02,"No large deletions",adj=0,cex=1.2,col=1)
	}
	if(MW[i]<0.1 & MW2[i]>=0.1 & LW[i]>=0.1 & MWMW2[i]>=0.1){
		text(4.7,m*0.02,"Potential deletion\nof OPN1MW exon5",adj=0,col=2,cex=1.2)
	}
	if(MW[i]>=0.1 & MW2[i]<0.1 & LW[i]>=0.1 & MWMW2[i]>=0.1){
		text(4.7,m*0.02,"Potential deletion\nof OPN1MW2 exon5",adj=0,col=2,cex=1.2)
	}

	nreads=sum(haplo[which(haplo[,1]==colnames(data3)[i]),3])

	text(4.7,m*0.56,"Exon 3 haplotypes:",adj=0,cex=1.2)
	x=haplo[which(haplo[,1]==colnames(data3)[i] & haplo[,3]>0),2:3]
	x=x[sort(as.numeric(x[,2]), index.return=TRUE, decreasing=TRUE)$ix,]
	for (j in 1:dim(x)[1]){
		text(4.7,m*(0.56-(j)*0.035),paste(x[j,1]," (n = ",x[j,2],")",sep=""),adj=0)
	}
	
	x=haploP[which(haploP[,1]==colnames(data3)[i] & as.numeric(haploP[,3])>1),2:3]
	np=length(which(haploP[,1]==colnames(data3)[i] & as.numeric(haploP[,3])>1))
	if(length(which(haploP[,1]==colnames(data3)[i] & as.numeric(haploP[,3])>1))>0){
		text(4.7,m*0.32,"Pathogenic haplotypes (>1 read):",adj=0,cex=1.2,col=2)
		if(np>1){
			x=x[sort(as.numeric(x[,2]), index.return=TRUE, decreasing=TRUE)$ix,]
			for (j in 1:dim(x)[1]){
				perc=round(as.numeric(x[j,2])/(nreads),digits=2)
				text(4.7,m*(0.32-(j)*0.035),paste(x[j,1]," n=",x[j,2],"/",nreads," (",perc*100,"%)",sep=""),adj=0,col=2)
				pa="Low"
				if(as.numeric(x[j,2])>1 & as.numeric(perc)>0.75){
				#if(as.numeric(x[j,1])>10 & as.numeric(perc)>0.75){
					pa="High"
					score1=score1+1
				}
				if(as.numeric(x[j,2])>1 & as.numeric(perc)<=0.75 & as.numeric(perc)>0.4){
				#if(as.numeric(x[j,1])>10 & as.numeric(perc)<=0.75 & as.numeric(perc)>0.4){
					score2=score2+0.5
				}
				if(as.numeric(x[j,2])>1 & as.numeric(perc)<=0.4 & as.numeric(perc)>0.1){
					score3=score3+0.5
				}
				events=rbind(events,c(colnames(data3)[i],sex[i],"Pathogenic haplotype",pa,x[j,1],x[j,2],nreads-as.numeric(x[j,2]),perc,"NA"))
			}
		}
		if(np==1){
			j=1
			perc=round(as.numeric(x[2])/(nreads),digits=2)
			text(4.7,m*(0.32-(j)*0.035),paste(x[j,1]," n=",x[j,2],"/",nreads," (",perc*100,"%)",sep=""),adj=0,col=2)
			pa="Low"
			if(as.numeric(x[2])>1 & as.numeric(perc)>0.75){
				pa="High"
				score1=score1+1
			}
			if(as.numeric(x[2])>1 & as.numeric(perc)<=0.75 & as.numeric(perc)>0.4){
			#if(as.numeric(x[j,1])>10 & as.numeric(perc)<=0.75 & as.numeric(perc)>0.4){
				score2=score2+0.5
			}
			if(as.numeric(x[2])>1 & as.numeric(perc)<=0.4 & as.numeric(perc)>0.1){
				score3=score3+0.5
			}
			events=rbind(events,c(colnames(data3)[i],sex[i],"Pathogenic haplotype",pa,x[j,1],x[j,2],nreads-as.numeric(x[2]),perc,"NA"))
		}
	} else {
		text(4.7,m*0.33,"No pathogenic\nhaplotypes detected",adj=0,cex=1.2,col=1)
	}

	NS=rareVariants[which(rareVariants[,1]==colnames(data3)[i]),]
	if(dim(NS)[1]>0){
		text(4.7,m*0.21,paste("Rare variant(s)*:",sep=""),adj=0,cex=1.2,col=2)
		for (j in 1:dim(NS)[1]){
			gnoMAX=0
			mut=as.numeric(NS[j,7])
			WT=as.numeric(NS[j,6])
			if(NS[j,11]!="NP_064445.2:p.?"){
				c=strsplit(NS[j,10],":")[[1]][2]
				gno=which(gnomAD[,7]==c)
				if(length(gno)==1){
					gno2=paste("gnomAD (",gnomAD[gno,1],"): ",signif(gnomAD[gno,10],digits=2),sep="")
					gnoMAX=max(gnoMAX,gnomAD[gno,10])
				}
				if(length(gno)>1){
					gno2=paste("gnomAD (",gnomAD[gno[1],1],"): ",signif(gnomAD[gno[1],10],digits=2),sep="")
					gnoMAX=max(gnoMAX,gnomAD[gno[1],10])
					for (k in 2:length(gno)){
						gno2=paste(gno2,", gnomAD (",gnomAD[gno[k],1],"): ",signif(gnomAD[gno[k],10],digits=2),sep="")
						gnoMAX=max(gnoMAX,gnomAD[gno[k],10])
					}
				}
				if(length(gno)==0){gno2="not in gnomAD"}

				t=paste("chrX:",NS[j,3],NS[j,4],">",NS[j,5],"\n",NS[j,11],sep="")
				t2=paste("chrX:",NS[j,3],NS[j,4],">",NS[j,5],", ",NS[j,10],", ",NS[j,11],sep="")
			} else {
				c=strsplit(NS[j,10],":")[[1]][2]
				gno=which(gnomAD[,7]==c)
				if(length(gno)==1){
					gno2=paste("gnomAD (",gnomAD[gno,1],"): ",signif(gnomAD[gno,10],digits=2),sep="")
					gnoMAX=max(gnoMAX,gnomAD[gno,10])
				}
				if(length(gno)>1){
					gno2=paste("gnomAD (",gnomAD[gno[1],1],"): ",signif(gnomAD[gno[1],10],digits=2),sep="")
					gnoMAX=max(gnoMAX,gnomAD[gno[1],10])
					for (k in 2:length(gno)){
						gno2=paste(gno2,", gnomAD (",gnomAD[gno[k],1],"): ",signif(gnomAD[gno[k],10],digits=2),sep="")
						gnoMAX=max(gnoMAX,gnomAD[gno[k],10])
					}
				}
				if(length(gno)==0){gno2="not in gnomAD"}
				t=paste("chrX:",NS[j,3],NS[j,4],">",NS[j,5],"\n",NS[j,10],sep="")
				t2=paste("chrX:",NS[j,3],NS[j,4],">",NS[j,5],", ",NS[j,10],sep="")
			}
			text(4.7,m*0.16-(j-1)*0.25,paste(t,"\n",mut,"/",WT+mut," reads (",round(100*mut/(WT+mut),digits=1),"%) ",sep=""),adj=0,cex=0.9,col=2)
			text(4.7,m*(-0.025),"* relative to OPN1LW gene but could affect other genes",adj=0,cex=0.45,col=1)
			pa="Low"
			#if(mut>10 & mut/(WT+mut)>0.75 & gnoMAX<0.0001){
			if(mut>1 & mut/(WT+mut)>0.75 & gnoMAX<0.0001){
				pa="High"
			}
			if(mut>1 & mut/(WT+mut)>0.75 & (NS[j,11]=="NP_064445.2:p.(Cys203Arg)" | NS[j,11]=="NP_064445.2:p.(Asn94Lys)" | NS[j,11]=="NP_064445.2:p.(Arg330Gln)" | NS[j,11]=="NP_064445.2:p.(Gly338Glu)" | NS[j,11]=="NP_064445.2:p.(Trp177Arg)" | grepl("Met1?",NS[j,11])==T  | grepl("Ter",NS[j,11])==T | grepl("+1G",NS[j,10])==T | grepl("+2T",NS[j,10])==T | grepl("-1G",NS[j,10])==T| grepl("-2A",NS[j,11])==T)){
			#if(mut>10 & mut/(WT+mut)>0.75 & (NS[j,11]=="NP_064445.2:p.(Cys203Arg)" | NS[j,11]=="NP_064445.2:p.(Asn94Lys)" | NS[j,11]=="NP_064445.2:p.(Arg330Gln)" | NS[j,11]=="NP_064445.2:p.(Gly338Glu)" | NS[j,11]=="NP_064445.2:p.(Trp177Arg)" | grepl("Met1?",NS[j,11])==T  | grepl("Ter",NS[j,11])==T | grepl("+1G",NS[j,10])==T | grepl("+2T",NS[j,10])==T | grepl("-1G",NS[j,10])==T| grepl("-2A",NS[j,11])==T)){
				score1=score1+1
			}
			if(mut>1 & mut/(WT+mut)>0.4 & mut/(WT+mut)<=0.75 & (NS[j,11]=="NP_064445.2:p.(Cys203Arg)" | NS[j,11]=="NP_064445.2:p.(Asn94Lys)" | NS[j,11]=="NP_064445.2:p.(Arg330Gln)" | NS[j,11]=="NP_064445.2:p.(Gly338Glu)" | NS[j,11]=="NP_064445.2:p.(Trp177Arg)" | grepl("Met1?",NS[j,11])==T  | grepl("Ter",NS[j,11])==T | grepl("+1G",NS[j,10])==T | grepl("+2T",NS[j,10])==T | grepl("-1G",NS[j,10])==T| grepl("-2A",NS[j,11])==T)){
			#if(mut>10 & mut/(WT+mut)>0.4 & mut/(WT+mut)<0.75 & (NS[j,11]=="NP_064445.2:p.(Cys203Arg)" | NS[j,11]=="NP_064445.2:p.(Asn94Lys)" | NS[j,11]=="NP_064445.2:p.(Arg330Gln)" | NS[j,11]=="NP_064445.2:p.(Gly338Glu)" | NS[j,11]=="NP_064445.2:p.(Trp177Arg)" | grepl("Met1?",NS[j,11])==T  | grepl("Ter",NS[j,11])==T | grepl("+1G",NS[j,10])==T | grepl("+2T",NS[j,10])==T | grepl("-1G",NS[j,10])==T| grepl("-2A",NS[j,11])==T)){
				score2=score2+1
			}
			if(mut>1 & mut/(WT+mut)>0.1 & mut/(WT+mut)<=0.4 & (NS[j,11]=="NP_064445.2:p.(Cys203Arg)" | NS[j,11]=="NP_064445.2:p.(Asn94Lys)" | NS[j,11]=="NP_064445.2:p.(Arg330Gln)" | NS[j,11]=="NP_064445.2:p.(Gly338Glu)" | NS[j,11]=="NP_064445.2:p.(Trp177Arg)" | grepl("Met1?",NS[j,11])==T  | grepl("Ter",NS[j,11])==T | grepl("+1G",NS[j,10])==T | grepl("+2T",NS[j,10])==T | grepl("-1G",NS[j,10])==T| grepl("-2A",NS[j,11])==T)){
				score3=score3+0.5
			}

			events=rbind(events,c(colnames(data3)[i],sex[i],"Rare variant",pa,t2,mut,WT,round(mut/(WT+mut),digits=2),gno2))
		}
	} else {
		text(4.7,m*0.175,paste("No coding variants\ndetected",sep=""),adj=0,cex=1.2,col=1)
	}

	if(min(exons[,i])<=0.1){
		sco=0
		
		# LCR
		if((medi2[1]>20 & exons[1,i]<0.1) | (medi2[1]>5 & exons[1,i]==0)){
			sco=max(sco,2)
		}
		
		for (k in 2:5){
			if((medi2[k]>20 & exons[k,i]<0.1) | (medi2[k]>5 & exons[k,i]==0)){
				sco=max(sco,0.5)
			}
		}
		score2=score2+sco

		sco=0
		for (k in 8:10){
			if((medi2[k]>20 & exons[k,i]<0.1) | (medi2[k]>5 & exons[k,i]==0)){
				sco=max(sco,0.5)
			}
		}
		if((medi2[6]+medi2[7]>20 & exons[6,i]+exons[7,i]<0.2) | (medi2[6]+medi2[7]>5 & exons[6,i]+exons[7,i]==0)){
			sco=max(sco,0.5)
		}
		score2=score2+sco

		
		sco=0
		for (k in 11:12){
			if((medi2[k]>20 & exons[k,i]<0.1) | (medi2[k]>5 & exons[k,i]==0)){
				sco=max(sco,1)
			}
		}
		score1=score1+sco
	}

	pheno="Normal"
	if(score3>=0.5){pheno="Normal / Color blindness"}
	if(score2==0.5){pheno="Color blindness"}
	if(score2>=1){pheno="Color blindness / severe"}
	if(score1>=1){pheno="Blue cone monochromacy / severe"}

	colp=3
	cexp=1.5
	if(pheno!="Normal"){colp=2; cexp=1.2}
	text(4.7,m*1,"Inferred phenotype: ",cex=1.5,adj=0)
	text(4.7,m*0.95,pheno,cex=cexp,adj=0,col=colp)

	# exon plot
	beg=12.8
	par(mgp=c(2.5,0.8,0),mar=c(8.1, 4.1, 4.1, 2.1))
	m=max(2,exons*1,na.rm=T)
	m=max(2.5,exons[,i]*1.2,na.rm=T)
	plot(1:4,ylim=c(0,m),xlim=c(1,15.5), xaxt = "n",xlab="",ylab="Coverage ratio vs. controls",cex=0,main=paste("Details on exons",sep=""),cex.main=2,cex.lab=1.4)
	
	for (j in 1:12){
		if(medi2[j]<10){
			rect(j-0.4,-0.05,j+0.4,m*1.02,border=NA,col="#ffcccc")
		}
		if(medi2[j]>=10 & medi2[j]<=20){
			rect(j-0.4,-0.05,j+0.4,m*1.02,border=NA,col="#fde7b5")
		}
	}

	for (j in 0:10){
		lines(c(-1,12.5),c(j,j),lty=2,col="darkgrey")
	}
	lines(c(-1,12.5),c(1,1))
	lines(c(-1,12.5),c(0,0))
	cols=rep(4,12)
	for (j in 1:12){
		if(exons[j,i]<0.1){cols[j]=2}
	}

	for (j in 1:12){
		rect(j-0.3,0,j+0.3,exons[j,i],col=cols[j])
		if(cols[j]!=4 & exons[j,i]<0.02){lines(c(j-0.3,j+0.3),c(exons[j,i],exons[j,i]),col=cols[j],lwd=6)}
	}

	n=dim(exons)[2]-1
	for (j in 1:12){
		points(rep(j,n)+runif(n,-0.2,0.2),exons[j,-i],pch=16,cex=0.5)
	}
	axis(1, at=1:12, labels=rownames(exons),las=2,cex.axis=1)

	text(beg,m*0.99,"Inferred sex:",cex=1.2,adj=0)
	text(beg,m*0.95,sex[i],cex=1,adj=0)
	legend(beg,m*0.91,c("No deletion","Deletion","Other individuals","Very low cov. (<10)","Low cov. (<20)"),col=c(4,2,1,"#ffcccc","#fde7b5"),pch=c(15,15,16,15,15),pt.cex=c(2,2,1,2,2),bg="white",y.intersp=1.5)

	if(min(exons[,i])>0.1){
		text(beg,m*0.6,"No deletions\ndetected",adj=0,col=1,cex=1.2)
	} else {
		text(beg,m*0.6,"Deletion(s) detected:",adj=0,col=2,cex=1.2)
		count=1
		for (j in c(1:5,8:12)){
			if(exons[j,i]<0.1){
				text(beg+0.2,m*(0.6-0.07*count),rownames(exons)[j],adj=0,col=2,cex=1)
				count=count+1
				pa="Low"
				if((medi2[j]>20 & exons[j,i]<0.1) | (medi2[j]>5 & exons[j,i]==0)){
					pa="High"
				}
				events=rbind(events,c(colnames(data3)[i],sex[i],"Hemi/homozygous deletion",pa,namesexons[j],"NA",exons2[j,i],round(exons[j,i],digits=2),"NA"))
			}
		}
		if(exons[6,i]<0.1 & exons[7,i]<0.1){
			pa="Low"
			if((medi2[6]+medi2[7]>20 & exons[6,i]+exons[7,i]<0.1) | (medi2[6]+medi2[7]>5 & exons[6,i]+exons[7,i]==0)){
				pa="High"
			}
			events=rbind(events,c(colnames(data3)[i],sex[i],"Hemi/homozygous deletion",pa,"OPN1MW/MW2-exon5","NA",exons2[6,i]+exons2[7,i],round(exons[6,i]+exons[7,i],digits=2),"NA"))
		}
		for (j in c(6,7)){
			if(exons[j,i]<0.1){
				text(beg+0.2,m*(0.6-0.07*count),rownames(exons)[j],adj=0,col=2,cex=1)
			}
		}

	}
	text(beg+0.2,m*0,"* LCR is usually not covered\n  and can have high variability.",cex=0.7,adj=0)

	phenotype=rbind(phenotype,c(colnames(data3)[i],sex[i],score1+score2/2,pheno))

	dev.off()
}

if(dim(events)[1]>0){
	for (i in 1:dim(events)[1]){
		events[i,3]=gsub("\n"," ",events[i,3])
	}
}

write.table(events,file=paste(dir,"/0.events.tsv",sep=""),quote=F,row.names=F,sep="\t")
write.table(phenotype[,c(1,2,4)],file=paste(dir,"/0.phenotypes-all.tsv",sep=""),quote=F,row.names=F,sep="\t")


