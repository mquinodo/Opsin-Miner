
#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

file1=args[1]
file2=args[2]
file3=args[3]
file4=args[4]
file5=args[5]
file6=args[6]
dir=args[7]
file8=args[8]
malesonly=args[9]

load(file6)

# reading files
data=read.table(file=file1,header=T,sep="\t",check.names = FALSE)
haplo=read.table(file=file2)
haploP=read.table(file=file3)
haploDNA=read.table(file=file4)
rareVariants=read.table(file=file5,header=T,row.names=NULL,sep="\t",colClasses = c("character"))
muts=read.table(file=file8)

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
MW=1:dim(data3)[2]
exons<-matrix(nrow=11,ncol=dim(data3)[2])
rownames(exons)=c("LCR*","ALL\nexon 3","ALL\nexon 6","OPN1LW\nexon 1","OPN1LW\nexon 2","OPN1LW\nexon 4","OPN1LW\nexon 5","OPN1MW\nexon 1","OPN1MW\nexon 2","OPN1MW\nexon 4","OPN1MW\nexon 5")
namesexons=c("LCR","ALL-exon3","ALLexon-6","OPN1LW-exon1","OPN1LW-exon2","OPN1LW-exon4","OPN1LW-exon5","OPN1MW-exon1","OPN1MW-exon2","OPN1MW-exon4","OPN1MW-exon5")
exons2=exons

# mean per target without sample with less than 10% of average
me=data3[,1]
for (i in 1:dim(data3)[1]){
	me[i]=mean(as.numeric(data3[i,which(as.numeric(data3[i,])>0.1*mean(as.numeric(data3[i,])))]))
}

data4=data[,-1]
for (i in 1:dim(data3)[2]){
	LW[i]=sum(data3[9:16,i])/sum(me[9:16])
	MW[i]=sum(data3[17:26,i])/sum(me[17:26])

	exons[1,i]=max(as.numeric(data3[3,i])/me[3]) # LCR
	exons[2,i]=max(as.numeric(data3[5,i])/me[5],as.numeric(data3[6,i])/me[6]) # exon 3
	exons[3,i]=max(as.numeric(data3[7,i])/me[7],as.numeric(data3[8,i])/me[8]) # exon 6
	exons[4,i]=max(as.numeric(data3[9,i])/me[9],as.numeric(data3[10,i])/me[10]) # LW exon1 
	exons[5,i]=max(as.numeric(data3[11,i])/me[11],as.numeric(data3[12,i])/me[12]) # LW exon2 
	exons[6,i]=max(as.numeric(data3[13,i])/me[13],as.numeric(data3[14,i])/me[14]) # LW exon4 
	exons[7,i]=max(as.numeric(data3[15,i])/me[15],as.numeric(data3[16,i])/me[16]) # LW exon5
	exons[8,i]=max(as.numeric(data3[21,i])/me[21],as.numeric(data3[22,i])/me[22]) # MW exon1 
	exons[9,i]=max(as.numeric(data3[23,i])/me[23],as.numeric(data3[24,i])/me[24]) # MW exon2 
	exons[10,i]=max(as.numeric(data3[25,i])/me[25],as.numeric(data3[26,i])/me[26]) # MW exon4 
	exons[11,i]=max(as.numeric(data3[17,i])/me[17],as.numeric(data3[18,i])/me[18])+max(as.numeric(data3[19,i])/me[19],as.numeric(data3[20,i])/me[20]) # MW exon5
	
	exons2[1,i]=max(as.numeric(data4[3,i])) # LCR
	exons2[2,i]=max(as.numeric(data4[5,i]),as.numeric(data4[6,i])) # exon 3
	exons2[3,i]=max(as.numeric(data4[7,i]),as.numeric(data4[8,i])) # exon 6
	exons2[4,i]=max(as.numeric(data4[9,i]),as.numeric(data4[10,i])) # LW exon1 
	exons2[5,i]=max(as.numeric(data4[11,i]),as.numeric(data4[12,i])) # LW exon2 
	exons2[6,i]=max(as.numeric(data4[13,i]),as.numeric(data4[14,i])) # LW exon4 
	exons2[7,i]=max(as.numeric(data4[15,i]),as.numeric(data4[16,i])) # LW exon5 
	exons2[8,i]=max(as.numeric(data4[21,i]),as.numeric(data4[22,i])) # MW exon1
	exons2[9,i]=max(as.numeric(data4[23,i]),as.numeric(data4[24,i])) # MW exon2
	exons2[10,i]=max(as.numeric(data4[25,i]),as.numeric(data4[26,i])) # MW exon4
	exons2[11,i]=max(as.numeric(data4[17,i]),as.numeric(data4[18,i]))+max(as.numeric(data4[19,i]),as.numeric(data4[20,i])) # MW exon5

}

medi1=c(median(apply(data4[9:16,],2,sum)),median(apply(data4[17:26,],2,sum)))
medi2=c(median(as.numeric(data4[3,])),median(apply(data4[5:6,],2,max)),median(apply(data4[7:8,],2,max)),
	median(apply(data4[9:10,],2,max)),median(apply(data4[11:12,],2,max)),median(apply(data4[13:14,],2,max)),median(apply(data4[15:16,],2,max)),
	median(apply(data4[21:22,],2,max)),median(apply(data4[23:24,],2,max)),median(apply(data4[25:26,],2,max)),median(apply(data4[17:18,],2,max)+apply(data4[19:20,],2,max)))

exons=exons[c(1,4:11,2:3),]
exons2=exons2[c(1,4:11,2:3),]
medi2=medi2[c(1,4:11,2:3)]
namesexons=namesexons[c(1,4:11,2:3)]

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
sex[which(z1<(-1.5))]="Male"
sex[which(z1>(1.5))]="Female"
sex[which(z1<=(1.5) & z1>=(-1.5))]="Undetermined"

pdf(file=paste(dir,"/07_plots/","0.sex.pdf",sep=""),height=10,width=20)
plot(z1)
dev.off()

if(malesonly=="Yes"){
	sex[1:length(sex)]="Male"
	sexout=cbind(colnames(data3),sex)
	colnames(sexout)=c("ID","Sex")
	write.table(sexout,file=paste(dir,"/0.sex.tsv",sep=""),quote=F,row.names=F,sep="\t")
} else {
	sexout=cbind(colnames(data3),sex,z1)
	colnames(sexout)=c("ID","Inferred-sex","Z-score")
	write.table(sexout,file=paste(dir,"/0.sex.tsv",sep=""),quote=F,row.names=F,sep="\t")
}

events<-matrix(ncol=9,nrow=0)
colnames(events)=c("#ID","Inferred-sex","Event","Pathogenicity","Details","Supporting-reads","WT-reads","Ratio","gnomAD-AF")

phenotype<-matrix(ncol=7,nrow=0)
colnames(phenotype)=c("#ID","Inferred-sex","Score","Phenotype","CN_OPN1LW","CN_OPN1MW","CN_hybrid")


# new analysis

d=data[,-1]
ALL<-matrix(ncol=14,nrow=dim(d)[2])
rownames(ALL)=colnames(data)[2:dim(data)[2]]
colnames(ALL)=c("Control-autosome","Control-chrX","LCR","LCR-control","Exon3","Exon6","LWExon1","LWExon2","LWExon4","LWExon5","MWExon1","MWExon2","MWExon4","MWExon5")
for (i in 1:dim(d)[2]){
    ALL[i,1]=max(d[1,i])/d[1,i]
    ALL[i,2]=max(d[2,i])/d[1,i]
    ALL[i,3]=max(d[3,i])/d[1,i]
    ALL[i,4]=max(d[4,i])/d[1,i]
    ALL[i,5]=max(d[5:6,i])/d[1,i]
    ALL[i,6]=max(d[7:8,i])/d[1,i]
    ALL[i,7]=max(d[9:10,i])/d[1,i]
    ALL[i,8]=max(d[11:12,i])/d[1,i]
    ALL[i,9]=max(d[13:14,i])/d[1,i]
    ALL[i,10]=max(d[15:16,i])/d[1,i]
    ALL[i,11]=max(d[21:22,i])/d[1,i]
    ALL[i,12]=max(d[23:24,i])/d[1,i]
    ALL[i,13]=max(d[25:26,i])/d[1,i]
    ALL[i,14]=(max(d[17:18,i])+max(d[19:20,i]))/d[1,i]
}

LWmin=as.numeric(apply(ALL[,7:10],1,min)*d[1,])
sel=which(LWmin>=10 & sex=="Male")

LWcor=cbind(ALL[,7]/mean(ALL[sel,7]),ALL[,8]/mean(ALL[sel,8]),ALL[,9]/mean(ALL[sel,9]),ALL[,10]/mean(ALL[sel,10]))
LWs=apply(LWcor,1,min)
LWs=LWs/mean(LWs[sel])

male=which(LWmin>=10 & sex=="Male")
female=which(LWmin>=10 & sex=="Female")
bad=which(LWmin<10)

MWmin=as.numeric(apply(ALL[,11:14],1,min))
sel=which(MWmin>=mean(MWmin)/10 & sex=="Male")

MWcor=cbind(ALL[,11]/mean(ALL[sel,11]),ALL[,12]/mean(ALL[sel,12]),ALL[,13]/mean(ALL[sel,13]),ALL[,14]/mean(ALL[sel,14]))
MWs=apply(MWcor,1,min)
MWs=MWs/mean(MWs[sel])

f=seq(median(MWs)*1,median(MWs)*2,0.01)
r=f
for (i in 1:length(f)){
    s1=abs(MWs*f[i]-1)
    s2=abs(MWs*f[i]-2)
    s3=abs(MWs*f[i]-3)
    s4=abs(MWs*f[i]-4)
    s5=abs(MWs*f[i]-5)
    s6=abs(MWs*f[i]-6)
    s7=abs(MWs*f[i]-7)
    s8=abs(MWs*f[i]-8)
    s9=abs(MWs*f[i]-9)
    s10=abs(MWs*f[i]-10)
    s=cbind(s1,s2,s3,s4,s5,s6,s7,s8,s9,s10)
    r[i]=sum(apply(s[sel,],1,min))
}
MWs=MWs*f[which.min(r)]

male=which(MWmin>=10 & sex=="Male")
female=which(MWmin>=10 & sex=="Female")
bad=which(MWmin<10)

cLW=round(LWs)
cMW=round(MWs)
sdLW=sd(LWs-cLW)
sdMW=sd(MWs-cMW)

LWmax=8
pLW<-matrix(nrow=length(cLW),ncol=LWmax)
for (i in 1:length(cLW)){
    tot=0:5
    tot[1]=dnorm(LWs[i],mean=0,sd=sdLW/4)
    for (j in 2:LWmax){
       tot[j]=dnorm(LWs[i],mean=j-1,sd=sdLW*j)
    }
    for (j in 1:LWmax){
       pLW[i,j]=tot[j]/sum(tot)
    }
}
cLW=apply(pLW,1,which.max)-1

pMW<-matrix(nrow=length(cMW),ncol=12)
for (i in 1:length(cMW)){
    tot=0:11
    tot[1]=dnorm(MWs[i],mean=0,sd=sdMW/4)
    for (j in 2:12){
       tot[j]=dnorm(MWs[i],mean=j-1,sd=sdMW*j)
    }
    for (j in 1:12){
       pMW[i,j]=tot[j]/sum(tot)
    }
}
cMW=apply(pMW,1,which.max)-1

#plot(apply(pLW,1,max),LWs)
#plot(apply(pMW,1,max),MWs)

out=cbind(cLW,cMW)
out2=apply(out, 1, paste, collapse="+")

male=which(sex=="Male")
female=which(sex=="Female")


for (i in 1:dim(pLW)[2]){
	pLW[,i]=round(pLW[,i],digits=2)
}

for (i in 1:dim(pMW)[2]){
	pMW[,i]=round(pMW[,i],digits=2)
}

colnames(pLW)=c("0 copy","1 copy","2 copies","3 copies","4 copies","5 copies","6 copies","7 copies")
colnames(pMW)=c("0 copy","1 copy","2 copies","3 copies","4 copies","5 copies","6 copies","7 copies","8 copies","9 copies","10 copies","11 copies")

x=as.matrix(rownames(ALL))
colnames(x)="ID"

write.table(cbind(x,pLW),file=paste(dir,"/06_copy-number/0.proba_CN-LW.tsv",sep=""),quote=F,row.names=F,sep="\t")
write.table(cbind(x,pMW),file=paste(dir,"/06_copy-number/0.proba_CN-MW.tsv",sep=""),quote=F,row.names=F,sep="\t")

cn=cbind(cLW,cMW)
colnames(cn)=c("CN_OPN1LW","CN_OPN1MW")

write.table(cbind(x,cn),file=paste(dir,"/0.copy-number.tsv",sep=""),quote=F,row.names=F,sep="\t")


ALL2=ALL

i=3 # LCR
LWmin=as.numeric(ALL[,i])
sel=which(LWmin>=mean(LWmin)/10 & sex=="Male")
ALL2[,i]=ALL[,i]/mean(ALL[sel,i])
ALL2[,i]=ALL2[,i]/mean(ALL2[sel,i])

for(i in 5:6){
	LWmin=as.numeric(ALL[,i])
	sel=which(LWmin>=mean(LWmin)/10 & sex=="Male")
	ALL2[,i]=ALL[,i]/mean(ALL[sel,i])
	ALL2[,i]=ALL2[,i]/mean(ALL2[sel,i])

	f=seq(median(ALL2[sel,i])*2,median(ALL2[sel,i])*3,0.01)
	r=f
	for (k in 1:length(f)){
	    s1=abs(ALL2[,i]*f[k]-1)
	    s2=abs(ALL2[,i]*f[k]-2)
	    s3=abs(ALL2[,i]*f[k]-3)
	    s4=abs(ALL2[,i]*f[k]-4)
	    s5=abs(ALL2[,i]*f[k]-5)
	    s6=abs(ALL2[,i]*f[k]-6)
	    s7=abs(ALL2[,i]*f[k]-7)
	    s8=abs(ALL2[,i]*f[k]-8)
	    s9=abs(ALL2[,i]*f[k]-9)
	    s10=abs(ALL2[,i]*f[k]-10)
	    s=cbind(s1,s2,s3,s4,s5,s6,s7,s8,s9,s10)
	    r[k]=sum(apply(s[sel,],1,min))
	}
	plot(f,r,ylim=c(0,8))
	ALL2[,i]=ALL2[,i]*f[which.min(r)]
}

for(i in 7:10){
	LWmin=as.numeric(ALL[,i])
	sel=which(LWmin>=mean(LWmin)/10 & sex=="Male")
	ALL2[,i]=ALL[,i]/mean(ALL[sel,i])
	ALL2[,i]=ALL2[,i]/mean(ALL2[sel,i])

	f=seq(median(ALL2[sel,i])*0.5,median(ALL2[sel,i])*1.5,0.01)
	r=f
	for (k in 1:length(f)){
	    s1=abs(ALL2[,i]*f[k]-1)
	    s2=abs(ALL2[,i]*f[k]-2)
	    s3=abs(ALL2[,i]*f[k]-3)
	    s4=abs(ALL2[,i]*f[k]-4)
	    s5=abs(ALL2[,i]*f[k]-5)
	    s6=abs(ALL2[,i]*f[k]-6)
	    s7=abs(ALL2[,i]*f[k]-7)
	    s8=abs(ALL2[,i]*f[k]-8)
	    s9=abs(ALL2[,i]*f[k]-9)
	    s10=abs(ALL2[,i]*f[k]-10)
	    s=cbind(s1,s2,s3,s4,s5,s6,s7,s8,s9,s10)
	    r[k]=sum(apply(s[sel,],1,min))
	}
	plot(f,r,ylim=c(0,8))
	ALL2[,i]=ALL2[,i]*f[which.min(r)]
}

for(i in 11:14){
	LWmin=as.numeric(ALL[,i])
	sel=which(LWmin>=mean(LWmin)/10 & sex=="Male")
	ALL2[,i]=ALL[,i]/mean(ALL[sel,i])
	ALL2[,i]=ALL2[,i]/mean(ALL2[sel,i])

	f=seq(median(ALL2[sel,i])*1,median(ALL2[sel,i])*2,0.01)
	r=f
	for (k in 1:length(f)){
	    s1=abs(ALL2[,i]*f[k]-1)
	    s2=abs(ALL2[,i]*f[k]-2)
	    s3=abs(ALL2[,i]*f[k]-3)
	    s4=abs(ALL2[,i]*f[k]-4)
	    s5=abs(ALL2[,i]*f[k]-5)
	    s6=abs(ALL2[,i]*f[k]-6)
	    s7=abs(ALL2[,i]*f[k]-7)
	    s8=abs(ALL2[,i]*f[k]-8)
	    s9=abs(ALL2[,i]*f[k]-9)
	    s10=abs(ALL2[,i]*f[k]-10)
	    s=cbind(s1,s2,s3,s4,s5,s6,s7,s8,s9,s10)
	    r[k]=sum(apply(s[sel,],1,min))
	}
	plot(f,r,ylim=c(0,8))
	ALL2[,i]=ALL2[,i]*f[which.min(r)]
}

ALL2=t(ALL2[,c(3,7,8,9,10,11,12,13,14,5,6)])
rownames(ALL2)=rownames(exons)


# plots
for (i in 1:dim(data3)[2]){

	score1=0
	score2=0
	score3=0
	score4=0
	score5=0
	score6=0

	n=length(LW)-1

	nreads=sum(haplo[which(haplo[,1]==colnames(data3)[i]),3])

	x=haplo[which(haplo[,1]==colnames(data3)[i] & haplo[,3]>1),2:3]
	x=x[sort(as.numeric(x[,2]), index.return=TRUE, decreasing=TRUE)$ix,]

	
	x=haploP[which(haploP[,1]==colnames(data3)[i] & as.numeric(haploP[,3])>1),2:3]
	np=length(which(haploP[,1]==colnames(data3)[i] & as.numeric(haploP[,3])>1))
	if(length(which(haploP[,1]==colnames(data3)[i] & as.numeric(haploP[,3])>1))>0){
		if(np>1){
			x=x[sort(as.numeric(x[,2]), index.return=TRUE, decreasing=TRUE)$ix,]
			for (j in 1:dim(x)[1]){
				perc=round(as.numeric(x[j,2])/(nreads),digits=2)
				pa="Low"
				if(as.numeric(x[j,2])>1 & as.numeric(perc)>0.75){
					pa="High"
					score1=score1+1
				}
				if(as.numeric(x[j,2])>1 & as.numeric(perc)<=0.75 & as.numeric(perc)>0.4){
					score2=score2+1
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
			pa="Low"
			if(as.numeric(x[2])>1 & as.numeric(perc)>0.75){
				pa="High"
				score1=score1+1
			}
			if(as.numeric(x[2])>1 & as.numeric(perc)<=0.75 & as.numeric(perc)>0.4){
				score2=score2+1
			}
			if(as.numeric(x[2])>1 & as.numeric(perc)<=0.4 & as.numeric(perc)>0.1){
				score3=score3+0.5
			}
			events=rbind(events,c(colnames(data3)[i],sex[i],"Pathogenic haplotype",pa,x[j,1],x[j,2],nreads-as.numeric(x[2]),perc,"NA"))
		}
	}

	NS=rareVariants[which(rareVariants[,1]==colnames(data3)[i]),]
	if(dim(NS)[1]>0){
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

				t=paste("chrX:",NS[j,3],NS[j,4],">",NS[j,5],"\n",NS[j,10],"\n",NS[j,11],sep="")
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
			pa="Low"
			pathoV=0
			if(mut>1 & mut/(WT+mut)>0.75 & gnoMAX<0.0001){
				pa="High"
			}
			if(mut>1 & mut/(WT+mut)>0.75 & (grepl("p.\\(Cys203Arg\\)",NS[j,11])==T | grepl("p.\\(Asn94Lys\\)",NS[j,11])==T | grepl("p.\\(Arg330Gln\\)",NS[j,11])==T | grepl("p.\\(Gly338Glu\\)",NS[j,11])==T | grepl("p.\\(Trp177Arg\\)",NS[j,11])==T | grepl("Met1\\?",NS[j,11])==T  | grepl("Ter",NS[j,11])==T | grepl("\\+1G",NS[j,10])==T | grepl("\\+2T",NS[j,10])==T | grepl("-1G",NS[j,10])==T| grepl("-2A",NS[j,11])==T)){
				score1=score1+1
				pathoV=1
			}
			if(mut>1 & grepl("p.\\(Trp177Arg\\)",NS[j,11])==T){
				score6=score6+1
			}
			if(mut>1 & mut/(WT+mut)>0.4 & mut/(WT+mut)<=0.75 & (grepl("p.\\(Cys203Arg\\)",NS[j,11])==T | grepl("p.\\(Asn94Lys\\)",NS[j,11])==T | grepl("p.\\(Arg330Gln\\)",NS[j,11])==T | grepl("p.\\(Gly338Glu\\)",NS[j,11])==T | grepl("p.\\(Trp177Arg\\)",NS[j,11])==T | grepl("Met1\\?",NS[j,11])==T  | grepl("Ter",NS[j,11])==T | grepl("\\+1G",NS[j,10])==T | grepl("\\+2T",NS[j,10])==T | grepl("-1G",NS[j,10])==T| grepl("-2A",NS[j,11])==T)){
				score2=score2+1
				pathoV=1
			}
			if(mut>1 & mut/(WT+mut)>0.1 & mut/(WT+mut)<=0.4 & (grepl("p.\\(Cys203Arg\\)",NS[j,11])==T | grepl("p.\\(Asn94Lys\\)",NS[j,11])==T | grepl("p.\\(Arg330Gln\\)",NS[j,11])==T | grepl("p.\\(Gly338Glu\\)",NS[j,11])==T | grepl("p.\\(Trp177Arg\\)",NS[j,11])==T | grepl("Met1\\?",NS[j,11])==T  | grepl("Ter",NS[j,11])==T | grepl("\\+1G",NS[j,10])==T | grepl("\\+2T",NS[j,10])==T | grepl("-1G",NS[j,10])==T| grepl("-2A",NS[j,11])==T)){
				score3=score3+0.5
				pathoV=1
			}
			if(pathoV==0){	
				events=rbind(events,c(colnames(data3)[i],sex[i],"Rare variant",pa,t2,mut,WT,round(mut/(WT+mut),digits=2),gno2))
				score4=score4+round(mut/(WT+mut),digits=2)
			}
			if(pathoV==1){	
				events=rbind(events,c(colnames(data3)[i],sex[i],"Pathogenic variant",pa,t2,mut,WT,round(mut/(WT+mut),digits=2),gno2))
			}
		}
	}

	if(min(exons[,i])<=0.1){
		sco=0
		
		# LCR
		if((medi2[1]>20 & exons[1,i]<0.1) | (medi2[1]>5 & exons[1,i]==0)){
			sco=max(sco,1)
		}
		score1=score1+sco
		
		# LW
		sco=0
		for (k in 2:5){
			if((medi2[k]>20 & exons[k,i]<0.1) | (medi2[k]>5 & exons[k,i]==0)){
				sco=max(sco,0.75)
			}
		}
		score5=score5+sco

		# MW
		sco=0
		for (k in 6:9){
			if((medi2[k]>20 & exons[k,i]<0.1) | (medi2[k]>5 & exons[k,i]==0)){
				sco=max(sco,0.75)
			}
		}
		if((medi2[6]+medi2[7]>20 & exons[6,i]+exons[7,i]<0.2) | (medi2[6]+medi2[7]>5 & exons[6,i]+exons[7,i]==0)){
			sco=max(sco,0.75)
		}
		score5=score5+sco

		# LW and MW
		sco=0
		for (k in 10:11){
			if((medi2[k]>20 & exons[k,i]<0.1) | (medi2[k]>5 & exons[k,i]==0)){
				sco=max(sco,1)
			}
		}
		score1=score1+sco
	}
	
	if(min(exons[,i])>0.1){
		
	} else {

		count=1
		for (j in c(1:11)){
			if(exons[j,i]<0.1){
				count=count+1
				pa="Low"
				if((medi2[j]>20 & exons[j,i]<0.1) | (medi2[j]>5 & exons[j,i]==0)){
					pa="High"
				}
				events=rbind(events,c(colnames(data3)[i],sex[i],"Hemi/homozygous deletion",pa,namesexons[j],"NA",exons2[j,i],round(exons[j,i],digits=2),"NA"))
			}
		}
	}

	pdf(file=paste(dir,"/07_plots/",rownames(ALL)[i],"_OPN1LW-MW_analysis.pdf",sep=""),height=7,width=14)
	par(mfrow=c(1,2))

	plot(1,xlim=c(0,10),ylim=c(0,9.5),xaxt ="n",yaxt ="n",xlab="",ylab="",cex=0,main=paste("Summary for ",colnames(ALL2)[i],sep=""),cex.main=1.5)
	text(0,9.3,"Graphical representation and copy-number (in total)",cex=1,adj=0)

	m=8.2
	lines(c(0,10),c(m,m))

	if((medi2[1]>20 & exons[1,i]<0.1) | (medi2[1]>5 & exons[1,i]==0)){
		rect(0.2,m-0.2,1,m+0.2,col="white",border="NA")
		lines(c(0.25,0.2,0.2,0.25),c(m-0.2,m-0.2,m+0.2,m+0.2))
		lines(c(0.95,1,1,0.95),c(m-0.2,m-0.2,m+0.2,m+0.2))
		text(0.6,m+0.5,"LCR deleted",col="red",cex=0.8)
	} else {
		rect(0.2,m-0.2,1,m+0.2,col="grey")
		text(0.6,m+0.5,"LCR",cex=0.8)
	}

	nLW=which.max(pLW[i,])
	nMW=which.max(pMW[i,])

	if(ALL[i,5]<mean(ALL[,5])/10 | ALL[i,6]<mean(ALL[,6])/10){
		nLW=1
		nMW=1
	}

	iLW=round(max(pLW[i,]),digits=2)
	iMW=round(max(pMW[i,]),digits=2)

	xl=7.3
	hybrid=0
	CH=0
	e1=0.292
	g1=0.15
	l1=0.3
	y1=m-l1/2
	y2=y1+l1
	x1=2

	fact=1.2

	poshyb=0
	if(ALL[i,5]>mean(ALL[,5])/10 & ALL[i,6]>mean(ALL[,6])/10 & (ALL[i,7]>mean(ALL[,7])/10 | ALL[i,11]>mean(ALL[,11])/10) & (ALL[i,8]>mean(ALL[,8])/10 | ALL[i,12]>mean(ALL[,12])/10) & (ALL[i,9]>mean(ALL[,9])/10 | ALL[i,13]>mean(ALL[,13])/10) & (ALL[i,10]>mean(ALL[,10])/10 | ALL[i,14]>mean(ALL[,14])/10) ){
		if((ALL[i,7]>mean(ALL[,7])/10 | ALL[i,8]>mean(ALL[,8])/10 | ALL[i,9]>mean(ALL[,9])/10 | ALL[i,10]>mean(ALL[,10])/10) & (ALL[i,11]>mean(ALL[,11])/10 | ALL[i,12]>mean(ALL[,12])/10 | ALL[i,13]>mean(ALL[,13])/10 | ALL[i,14]>mean(ALL[,14])/10)){
			poshyb=1
		}
	}

	if(nLW==1 & nMW>1 & poshyb==1){
		colo=c("2","2","2","2","2","2")
		nMWtest=min(which(pMW[i,]>0.2))
		if(ALL[i,7]<mean(ALL[,7])/10 & ALL2[6,i]>fact*(nMWtest-1)){colo[1]=3}
		if(ALL[i,8]<mean(ALL[,8])/10 & ALL2[7,i]>fact*(nMWtest-1)){colo[2]=3}
		if(ALL[i,9]<mean(ALL[,9])/10 & ALL2[8,i]>fact*(nMWtest-1)){colo[4]=3}
		if(ALL[i,10]<mean(ALL[,10])/10 & ALL2[9,i]>fact*(nMWtest-1)){colo[5]=3; colo[6]=3}
		if(colo[2]!=colo[4]){colo[3]="yellow"}
		if(colo[2]==3 & colo[4]==3){colo[3]=3}
		if(colo[1]!="2" | colo[2]!="2" | colo[3]!="2" | colo[4]!="2" | colo[5]!="2" | colo[6]!="2"){
			hybrid=1
			for(d in 1:6){
				rect(x1+e1*(d-1)+g1*(d-1),y1,x1+e1*d+g1*(d-1),y2,col=colo[d])
				text(x1+e1*(d-1)+g1*(d-1)+e1/2,m,d,cex=0.8)
			}

			ok=c()
			if(ALL[i,7]>mean(ALL[,7])/10){ok=c(ok,2)}
			if(ALL[i,8]>mean(ALL[,8])/10){ok=c(ok,3)}
			if(ALL[i,9]>mean(ALL[,9])/10){ok=c(ok,4)}
			if(ALL[i,10]>mean(ALL[,10])/10){ok=c(ok,5)}
			CH=round(mean(ALL2[ok,i]))

			colo2=colo
			colo2[which(colo2==2)]="Red"
			colo2[which(colo2==3)]="Green"
			colo2[which(colo2=="yellow")]="NA"
			colo2=paste(colo2,collapse=",")
			events=rbind(events,c(colnames(data3)[i],sex[i],"Hybrid gene","High",colo2,"NA","NA","NA","NA"))
		}
	}

	if(nLW>1 & nMW==1 & poshyb==1){
		colo=c(3,3,3,3,3,3)
		nLWtest=min(which(pLW[i,]>0.2))
		if(ALL[i,11]<mean(ALL[,11])/10 & ALL2[2,i]>fact*(nLWtest-1)){colo[1]=2}
		if(ALL[i,12]<mean(ALL[,12])/10 & ALL2[3,i]>fact*(nLWtest-1)){colo[2]=2}
		if(ALL[i,13]<mean(ALL[,13])/10 & ALL2[4,i]>fact*(nLWtest-1)){colo[4]=2}
		if(ALL[i,14]<mean(ALL[,14])/10 & ALL2[5,i]>fact*(nLWtest-1)){colo[5]=2; colo[6]=2}
		if(colo[2]!=colo[4]){colo[3]="yellow"}
		if(colo[2]==2 & colo[4]==2){colo[3]=2}
		if(colo[1]!="3" | colo[2]!="3" | colo[3]!="3" | colo[4]!="3" | colo[5]!="3" | colo[6]!="3"){
			hybrid=1
			x1=6
			for(d in 1:6){
				rect(x1+e1*(d-1)+g1*(d-1),y1,x1+e1*d+g1*(d-1),y2,col=colo[d])
				text(x1+e1*(d-1)+g1*(d-1)+e1/2,m,d,cex=0.8)
			}

			ok=c()
			if(ALL[i,11]>mean(ALL[,11])/10){ok=c(ok,6)}
			if(ALL[i,12]>mean(ALL[,12])/10){ok=c(ok,7)}
			if(ALL[i,13]>mean(ALL[,13])/10){ok=c(ok,8)}
			if(ALL[i,14]>mean(ALL[,14])/10){ok=c(ok,9)}
			CH=round(mean(ALL2[ok,i]))

			colo2=colo
			colo2[which(colo2==2)]="Red"
			colo2[which(colo2==3)]="Green"
			colo2[which(colo2=="yellow")]="NA"
			colo2=paste(colo2,collapse=",")
			events=rbind(events,c(colnames(data3)[i],sex[i],"Hybrid gene","High",colo2,"NA","NA","NA","NA"))
		}
	}

	if(nLW==1 & nMW==1 & poshyb==1){
		colo=c(2,2,2,2,2,2)
		if(ALL[i,7]<mean(ALL[,7])/10){colo[1]=3}
		if(ALL[i,8]<mean(ALL[,8])/10){colo[2]=3}
		if(ALL[i,9]<mean(ALL[,9])/10){colo[4]=3}
		if(ALL[i,10]<mean(ALL[,10])/10){colo[5]=3; colo[6]=3}
		if(colo[2]!=colo[4]){colo[3]="yellow"}
		if(colo[2]==3 & colo[4]==3){colo[3]=3}
		hybrid=1
		CH=1
		for(d in 1:6){
			rect(x1+e1*(d-1)+g1*(d-1),y1,x1+e1*d+g1*(d-1),y2,col=colo[d])
			text(x1+e1*(d-1)+g1*(d-1)+e1/2,m,d,cex=0.8)
		}

		colo2=colo
		colo2[which(colo2==2)]="Red"
		colo2[which(colo2==3)]="Green"
		colo2[which(colo2=="yellow")]="NA"
		colo2=paste(colo2,collapse=",")
		events=rbind(events,c(colnames(data3)[i],sex[i],"Hybrid gene","High",colo2,"NA","NA","NA","NA"))
	}


	x1=2

	if(nLW>1){
		x1=2
		x2=4
		x3=4.5
		l1=0.3
		l2=0.2
		y1=m-l1/2
		y2=y1+l1
		y3=y1-l2
		y4=y2+l2
		y5=(y1+y2)/2
		polygon(c(x1,x1,x2,x2,x3,x2,x2),c(y1,y2,y2,y4,y5,y3,y1),col=2,border="NA")
		lines(c(x1,x1,x2,x2,x3,x2,x2,x1),c(y1,y2,y2,y4,y5,y3,y1,y1))
		text(x1+(x3-x1)/2,m+0.5,"OPN1LW",cex=0.8)
		if(nLW==2){
			text(x1+(x3-x1)/2,m-0.5,paste(nLW-1," copy (p=",pLW[i,nLW],")",sep=""),cex=0.8)
		}
		if(nLW<=4 & nLW>2){
			text(x1+(x3-x1)/2,m-0.5,paste(nLW-1," copies (p=",pLW[i,nLW],")",sep=""),cex=0.8)
		}
		if(nLW>4){
			text(x1+(x3-x1)/2,m-0.5,paste(">3 copies",sum(pLW[i,nLW:dim(pLW)[2]]),")",sep=""),cex=0.8)
		}
	}

	if(nLW==1 & hybrid==0){

		colo=c(2,2,2,2,2,2)
		nLWtest=min(which(pLW[i,]>0.2))
		nMWtest=min(which(pMW[i,]>0.2))

		good=c(1:6)
		bad=c()
		if(ALL[i,7]<mean(ALL[,7])/10){bad=c(bad,1); good=good[-which(good==1)]}
		if(ALL[i,8]<mean(ALL[,8])/10){bad=c(bad,2); good=good[-which(good==2)]}
		if(ALL[i,9]<mean(ALL[,9])/10){bad=c(bad,4); good=good[-which(good==4)]}
		if(ALL[i,10]<mean(ALL[,10])/10){bad=c(bad,5); good=good[-which(good==5)]}
		if(ALL[i,10]<mean(ALL[,10])/10 & ALL2[11,i]-1<fact*(nLWtest+nMWtest-2)){bad=c(bad,6); good=good[-which(good==6)];}
		if(ALL[i,6]<mean(ALL[,6])/10){bad=c(bad,6); good=good[which(good!=6)];}
		if(ALL[i,8]<mean(ALL[,8])/10 | ALL[i,9]<mean(ALL[,9])/10){bad=c(bad,3); good=good[-which(good==3)]}

		for(d in bad){
			rect(x1+e1*(d-1)+g1*(d-1),y1,x1+e1*d+g1*(d-1),y2,col="white",border="NA")
			lines(c(x1+e1*(d-1)+g1*(d-1)+0.05,x1+e1*(d-1)+g1*(d-1),x1+e1*(d-1)+g1*(d-1),x1+e1*(d-1)+g1*(d-1)+0.05),c(m-0.2,m-0.2,m+0.2,m+0.2))
			lines(c(x1+e1*d+g1*(d-1)-0.05,x1+e1*d+g1*(d-1),x1+e1*d+g1*(d-1),x1+e1*d+g1*(d-1)-0.05),c(m-0.2,m-0.2,m+0.2,m+0.2))
		}
		for(d in good){
			rect(x1+e1*(d-1)+g1*(d-1),y1,x1+e1*d+g1*(d-1),y2,col=colo[d])
			text(x1+e1*(d-1)+g1*(d-1)+e1/2,m,d,cex=0.8)
		}

		# whole del
		if(ALL[i,7]<mean(ALL[,7])/10 & ALL[i,8]<mean(ALL[,8])/10 & ALL[i,9]<mean(ALL[,9])/10 & ALL[i,10]<mean(ALL[,10])/10){
			rect(2,m-0.22,4.5,m+0.22,col="white",border="NA")
			lines(c(4.45,4.5,4.5,4.45),c(m-0.2,m-0.2,m+0.2,m+0.2))
			lines(c(2.05,2,2,2.05),c(m-0.2,m-0.2,m+0.2,m+0.2))
			text(3.25,m+0.5,"OPN1LW deleted",col="red",cex=0.8)
		}
		if(length(good)>0){
			text(3.25,m+0.5,"Partial OPN1LW",col="red",cex=0.8)
		}

	}

	if(nLW==1 & hybrid==1){
		x1=2
		x2=4
		x3=4.5
		l1=0.3
		l2=0.2
		y1=m-l1/2
		y2=y1+l1
		y3=y1-l2
		y4=y2+l2
		y5=(y1+y2)/2
		text(x1+(x3-x1)/2,m+0.5,"Hybrid gene",cex=0.8)
		if(CH>1){
			text(x1+(x3-x1)/2,m-0.5,paste(CH," copies",sep=""),cex=0.8)
		}
		if(CH==1){
			text(x1+(x3-x1)/2,m-0.5,paste(CH," copy",sep=""),cex=0.8)
		}
	}

	if(nMW==1 & hybrid==0){

		x1=6

		colo=c(3,3,3,3,3,3)
		nLWtest=min(which(pLW[i,]>0.2))
		nMWtest=min(which(pMW[i,]>0.2))

		good=c(1:6)
		bad=c()
		if(ALL[i,11]<mean(ALL[,11])/10){bad=c(bad,1); good=good[-which(good==1)]}
		if(ALL[i,12]<mean(ALL[,12])/10){bad=c(bad,2); good=good[-which(good==2)]}
		if(ALL[i,13]<mean(ALL[,13])/10){bad=c(bad,4); good=good[-which(good==4)]}
		if(ALL[i,14]<mean(ALL[,14])/10){bad=c(bad,5); good=good[-which(good==5)]}
		if(ALL[i,14]<mean(ALL[,14])/10 & ALL2[11,i]-1<fact*(nLWtest+nMWtest-2)){bad=c(bad,6); good=good[-which(good==6)];}
		if(ALL[i,6]<mean(ALL[,6])/10){bad=c(bad,6); good=good[which(good!=6)];}
		if(ALL[i,12]<mean(ALL[,12])/10 | ALL[i,13]<mean(ALL[,13])/10){bad=c(bad,3); good=good[-which(good==3)]}

		for(d in bad){
			rect(x1+e1*(d-1)+g1*(d-1),y1,x1+e1*d+g1*(d-1),y2,col="white",border="NA")
			lines(c(x1+e1*(d-1)+g1*(d-1)+0.05,x1+e1*(d-1)+g1*(d-1),x1+e1*(d-1)+g1*(d-1),x1+e1*(d-1)+g1*(d-1)+0.05),c(m-0.2,m-0.2,m+0.2,m+0.2))
			lines(c(x1+e1*d+g1*(d-1)-0.05,x1+e1*d+g1*(d-1),x1+e1*d+g1*(d-1),x1+e1*d+g1*(d-1)-0.05),c(m-0.2,m-0.2,m+0.2,m+0.2))
		}
		for(d in good){
			rect(x1+e1*(d-1)+g1*(d-1),y1,x1+e1*d+g1*(d-1),y2,col=colo[d])
			text(x1+e1*(d-1)+g1*(d-1)+e1/2,m,d,cex=0.8)
		}

		# whole del
		if(ALL[i,11]<mean(ALL[,11])/10 & ALL[i,12]<mean(ALL[,12])/10 & ALL[i,13]<mean(ALL[,13])/10 & ALL[i,14]<mean(ALL[,14])/10){
			rect(6,m-0.22,8.5,m+0.22,col="white",border="NA")
			lines(c(8.45,8.5,8.5,8.45),c(m-0.2,m-0.2,m+0.2,m+0.2))
			lines(c(6.05,6,6,6.05),c(m-0.2,m-0.2,m+0.2,m+0.2))
			text(7.25,m+0.5,"OPN1MW deleted",col="red",cex=0.8)
		}
		if(length(good)>0){
			text(7.25,m+0.5,"Partial OPN1MW",col="red",cex=0.8)
		}

	}

	if(nMW==1 & hybrid==2){
		x1=6
		x2=8
		x3=8.5
		l1=0.3
		l2=0.2
		y1=m-l1/2
		y2=y1+l1
		y3=y1-l2
		y4=y2+l2
		y5=(y1+y2)/2
		# polygon(c(x1,x1,x2,x2,x3,x2,x2),c(y1,y2,y2,y4,y5,y3,y1),col=3,border="NA")
		# polygon(c(x1,x1,x1+(x3-x1)/3*2,x1+(x3-x1)/3),c(y1,y2,y2,y1),col=2,border="NA")
		# lines(c(x1,x1,x2,x2,x3,x2,x2,x1),c(y1,y2,y2,y4,y5,y3,y1,y1))
		text(x1+(x3-x1)/2,m+0.5,"Hybrid gene",cex=0.8)
		if(CH>1){
			text(x1+(x3-x1)/2,m-0.5,paste(CH," copies",sep=""),cex=0.8)
		}
		if(CH==1){
			text(x1+(x3-x1)/2,m-0.5,paste(CH," copy",sep=""),cex=0.8)
		}
	}

	if(nMW>1){
		x1=6
		x2=8
		x3=8.5
		l1=0.3
		l2=0.2
		y1=m-l1/2
		y2=y1+l1
		y3=y1-l2
		y4=y2+l2
		y5=(y1+y2)/2
		polygon(c(x1,x1,x2,x2,x3,x2,x2),c(y1,y2,y2,y4,y5,y3,y1),col=3,border="NA")
		lines(c(x1,x1,x2,x2,x3,x2,x2,x1),c(y1,y2,y2,y4,y5,y3,y1,y1))
		text(x1+(x3-x1)/2,m+0.5,"OPN1MW",cex=0.8)
		if(nMW==2){
			text(x1+(x3-x1)/2,m-0.5,paste(nMW-1," copy (p=",pMW[i,nMW],")",sep=""),cex=0.8)
		}
		if(nMW<=6 & nMW>2){
			text(x1+(x3-x1)/2,m-0.5,paste(nMW-1," copies (p=",pMW[i,nMW],")",sep=""),cex=0.8)
		}
		if(nMW>6){
			text(x1+(x3-x1)/2,m-0.5,paste(">5 copies (p=",sum(pMW[i,nMW:dim(pMW)[2]]),")",sep=""),cex=0.8)
		}
	}

	xl=7
	xlt=7

	nreads=sum(haplo[which(haplo[,1]==colnames(data3)[i]),3])
	text(0,xl,"Exon 3 protein haplotypes:",adj=0,cex=1)
	x=haplo[which(haplo[,1]==colnames(data3)[i]),2:3]
	x=x[sort(as.numeric(x[,2]), index.return=TRUE, decreasing=TRUE)$ix,]
	if(dim(x)[1]>0){
		for (j in 1:dim(x)[1]){
			if(x[j,1]=="LIAVA" | x[j,1]=="LVAVA" | x[j,1]=="MIAVA" | x[j,1]=="LIAVS" | x[j,1]=="LIVVA"){
					perc=round(as.numeric(x[j,2])/(nreads),digits=2)
					text(0,xl-0.4,paste(x[j,1]," n=",x[j,2],"/",nreads," (",perc*100,"%, pathogenic)",sep=""),adj=0,cex=0.8,col="red")
				} else {
					text(0,xl-0.4,paste(x[j,1]," (n = ",x[j,2],")",sep=""),adj=0,cex=0.8)
				}		
			xl=xl-0.4
		}
	}
	xl=xl+0.15

	nreads=sum(haploDNA[which(haploDNA[,1]==colnames(data3)[i]),3])
	text(5,xlt,"Exon 3 DNA haplotypes:",adj=0,cex=1)
	x=haploDNA[which(haploDNA[,1]==colnames(data3)[i]),2:3]
	x=x[sort(as.numeric(x[,2]), index.return=TRUE, decreasing=TRUE)$ix,]
	for (j in 1:dim(x)[1]){
		text(5,xlt-0.4,paste(x[j,1]," (n = ",x[j,2],")",sep=""),adj=0,cex=0.8)
		xlt=xlt-0.4
	}
	xlt=xlt+0.15

	xl=min(xl,xlt)-1
	

	NS=rareVariants[which(rareVariants[,1]==colnames(data3)[i]),]
	if(dim(NS)[1]>0){
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

				t1=paste("chrX:",NS[j,3],NS[j,4],">",NS[j,5],sep="")
				t2=paste(NS[j,10],", ",NS[j,11],sep="")
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
				t1=paste("chrX:",NS[j,3],NS[j,4],">",NS[j,5],sep="")
				t2=paste(NS[j,10],sep="")
			}
			pathoV=0
			if(mut>1 & mut/(WT+mut)>0.75 & ( grepl("p.\\(Cys203Arg\\)",NS[j,11])==T | grepl("p.\\(Asn94Lys\\)",NS[j,11])==T | grepl("p.\\(Arg330Gln\\)",NS[j,11])==T | grepl("p.\\(Gly338Glu\\)",NS[j,11])==T | grepl("p.\\(Trp177Arg\\)",NS[j,11])==T | grepl("Met1\\?",NS[j,11])==T  | grepl("Ter",NS[j,11])==T | grepl("\\+1G",NS[j,10])==T | grepl("\\+2T",NS[j,10])==T | grepl("-1G",NS[j,10])==T| grepl("-2A",NS[j,11])==T)){
				pathoV=1
			}
			if(mut>1 & mut/(WT+mut)>0.4 & mut/(WT+mut)<=0.75 & (grepl("p.\\(Cys203Arg\\)",NS[j,11])==T | grepl("p.\\(Asn94Lys\\)",NS[j,11])==T | grepl("p.\\(Arg330Gln\\)",NS[j,11])==T | grepl("p.\\(Gly338Glu\\)",NS[j,11])==T | grepl("p.\\(Trp177Arg\\)",NS[j,11])==T  | grepl("Met1\\?",NS[j,11])==T  | grepl("Ter",NS[j,11])==T | grepl("\\+1G",NS[j,10])==T | grepl("\\+2T",NS[j,10])==T | grepl("-1G",NS[j,10])==T| grepl("-2A",NS[j,11])==T)){
				pathoV=1
			}
			if(mut>1 & mut/(WT+mut)>0.1 & mut/(WT+mut)<=0.4 & (grepl("p.\\(Cys203Arg\\)",NS[j,11])==T | grepl("p.\\(Asn94Lys\\)",NS[j,11])==T | grepl("p.\\(Arg330Gln\\)",NS[j,11])==T | grepl("p.\\(Gly338Glu\\)",NS[j,11])==T | grepl("p.\\(Trp177Arg\\)",NS[j,11])==T | grepl("Met1\\?",NS[j,11])==T  | grepl("Ter",NS[j,11])==T | grepl("\\+1G",NS[j,10])==T | grepl("\\+2T",NS[j,10])==T | grepl("-1G",NS[j,10])==T| grepl("-2A",NS[j,11])==T)){
				pathoV=1
			}
			if(pathoV==0){	
				text(0,xl,paste("Rare variant (AF<1E-5, VUS, hg19)*:",sep=""),adj=0,cex=1,col="orange")
				text(0,xl-0.3,paste(t1,", ",mut,"/",WT+mut," reads (",round(100*mut/(WT+mut),digits=1),"%) ","\n",t2,sep=""),adj=c(0,1),cex=0.8,col="orange")
				#text(0,xl-1.2,"* relative to OPN1LW gene but could affect other genes",adj=0,cex=0.6,col=1)
				if(NS[j,10]=="NM_020061.5:c.619G>A" | NS[j,10]=="NM_000513.2:c.619G>A"){
					Lmut=muts[which(muts[,1]==colnames(data3)[i] & muts[,2]=="V207M-LW"),3]
					Mmut=muts[which(muts[,1]==colnames(data3)[i] & muts[,2]=="V207M-MW"),3]
					Lwt=muts[which(muts[,1]==colnames(data3)[i] & muts[,2]=="V207-LW"),3]
					Mwt=muts[which(muts[,1]==colnames(data3)[i] & muts[,2]=="V207-MW"),3]
					text(0,xl-0.9,  paste(Lmut,"/",Lwt+Lmut," reads (",round(100*Lmut/(Lmut+Lwt+0.00001),digits=1),"%) on OPN1LW-exon4 and ",Mmut,"/",Mwt+Mmut," reads (",round(100*Mmut/(Mwt+Mmut+0.00001),digits=1),"%) on OPN1MW-exon4",sep="")  ,adj=c(0,1),cex=0.8,col="orange")
					xl=xl-0.3
				}
				xl=xl-1.4
			}
			if(pathoV==1){	
				text(0,xl,paste("Pathogenic variant (hg19)*:",sep=""),adj=0,cex=1,col="red")
				text(0,xl-0.3,paste(t1,", ",mut,"/",WT+mut," reads (",round(100*mut/(WT+mut),digits=1),"%) ","\n",t2,sep=""),adj=c(0,1),cex=0.8,col="red")
				#text(0,xl-1.2,"* relative to OPN1LW gene but could affect other genes",adj=0,cex=0.6,col=1)
				if(NS[j,10]=="NM_020061.5:c.607T>C" | NS[j,10]=="NM_000513.2:c.607T>C"){
					Lmut=muts[which(muts[,1]==colnames(data3)[i] & muts[,2]=="C203R-LW"),3]
					Mmut=muts[which(muts[,1]==colnames(data3)[i] & muts[,2]=="C203R-MW"),3]
					Lwt=muts[which(muts[,1]==colnames(data3)[i] & muts[,2]=="C203-LW"),3]
					Mwt=muts[which(muts[,1]==colnames(data3)[i] & muts[,2]=="C203-MW"),3]
					text(0,xl-0.9,  paste(Lmut,"/",Lwt+Lmut," reads (",round(100*Lmut/(Lmut+Lwt+0.00001),digits=1),"%) on OPN1LW-exon4 and ",Mmut,"/",Mwt+Mmut," reads (",round(100*Mmut/(Mwt+Mmut+0.00001),digits=1),"%) on OPN1MW-exon4",sep="")  ,adj=c(0,1),cex=0.8,col="red")
					xl=xl-0.3
				}
				
				xl=xl-1.4
			}
		}
	} else {
		text(0,xl,paste("No coding variants detected",sep=""),adj=0,cex=1,col=1)
		xl=xl-0.8
	}



	if(min(exons[,i])>0.1){
		text(0,xl,"No deletions detected",adj=0,col=1)
	} else {
		text(0,xl,"Deletion(s) detected:",adj=0,col="red",cex=1)
		xl=xl-0
		count=1
		for (j in c(1:11)){
			if(exons[j,i]<0.1){
				text(0,xl-0.4,namesexons[j],adj=0,col="red",cex=0.8)
				xl=xl-0.4
			}
		}
	}
	#text(0,0,"# LCR is not covered by WES and can have high variability.",cex=0.7,adj=0)

	# score6 Trp177Arg
	# score5 dels L or M
	# score4 = rare VUS
	# score3 = events below 40% (0.5 each)
	# score2 = events 40-75% (1 each)
	# score1 = events above 75% (1 each)

	pheno="Normal"
	if(score4>=0.1){pheno="Normal with VUS"} 
	if(score4>=0.4){pheno="Normal / Color vision deficiency possible (VUS)"} 
	if(score4>=0.75){pheno="Normal / Color vision deficiency possible (VUS)"}

	if(score3>=0.5){pheno="Inconclusive"}
	if(score2>=1){pheno="Inconclusive"}
	if(score5>=0.75){pheno="Color vision deficiency suggested"}
	if(score5>=1.5 & hybrid==0){pheno="BCM suggested"}
	if(score5>=1.5 & hybrid>0){pheno="Color vision deficiency suggested"}
	if(score1>=1){pheno="BCM suggested"} 

	if(score3>=0.5 & score4>=0.1){pheno="Inconclusive (with VUS)"}
	if(score2>=1 & score4>=0.1){pheno="Inconclusive (with VUS)"}
	if(score5>=0.75 & score4>=0.1){pheno="Color vision deficiency suggested (VUS)"}
	if(score5>=1.5 & score4>=0.1){pheno="Color vision deficiency suggested (VUS)"}
	if(score5>=1.5 & score1>=1){pheno="BCM suggested (VUS)"}
	if(score5>=1.5 & score4>=0.75){pheno="Color vision deficiency or BCM (VUS)"}
	if(score1>=1 & score4>=0.1){pheno="BCM suggested (VUS)"}

	if(score6>=1){pheno="Cone dystrophy suggested"}

	phenotype=rbind(phenotype,c(colnames(data3)[i],sex[i],score1+score2/2,pheno,nLW-1,nMW-1,CH))


	p=pheno

	if(p=="Normal" | p=="Normal with VUS"){
		text(5,xl-1,paste("Inferred phenotype: ",p,sep=""),col=3,cex=1.1)
	}
	if(p!="Normal" & p!="Normal with VUS"){
		text(5,xl-1,paste("Inferred phenotype: ",p,sep=""),col="red",cex=1.1)
	}

	if(abs( data[1,i+1]-mean(as.numeric(data[1,2:dim(data)[2]]))) / sd(as.numeric(data[1,2:dim(data)[2]])) >=3 ){
		text(5,-1,paste("Low quality, results not reliable!",sep=""),col="red",cex=1.2, xpd=NA)
	}
	if(sex[i]=="Female"){
		text(0,-0.55,paste("Warning: Hybrid gene prediction and copy-numbers in females is not fully reliable.",sep=""),adj=0,col="red",cex=0.7, xpd=NA)
	}
	if(sex[i]=="Undetermined"){
		text(0,-0.55,paste("Warning: Sex could not be predicted. Results can be unreliable.",sep=""),col="red",cex=0.7, xpd=NA,adj=0)
	}


	# exon plot
	beg=11.8
	#par(mgp=c(2.5,0.8,0),mar=c(8.1, 4.1, 4.1, 2.1))
	m=max(2,ALL2*1,na.rm=T)
	m=max(2.5,ALL2[,i]*1.2,na.rm=T)
	m=max(m,5)
	plot(1:4,ylim=c(0,m),xlim=c(1,14.5), xaxt = "n",xlab="",ylab="Inferred copy-number",cex=0,main=paste("Copy-number of exons and deletions",sep=""),cex.main=1.5,cex.lab=1)
	
	for (j in 1:11){
		if(medi2[j]<10){
			rect(j-0.4,-0.05,j+0.4,m*1.02,border=NA,col="#ffcccc")
		}
		if(medi2[j]>=10 & medi2[j]<=20){
			rect(j-0.4,-0.05,j+0.4,m*1.02,border=NA,col="#fde7b5")
		}
	}

	for (j in 0:20){
		lines(c(-1,11.5),c(j,j),lty=2,col="darkgrey")
	}
	lines(c(-1,11.5),c(0,0))
	cols=rep(4,11)
	for (j in 1:11){
		if(ALL2[j,i]<0.1){cols[j]=2}
	}

	for (j in 1:11){
		rect(j-0.3,0,j+0.3,ALL2[j,i],col=cols[j])
		if(cols[j]!=4 & ALL2[j,i]<0.02){lines(c(j-0.3,j+0.3),c(ALL2[j,i],ALL2[j,i]),col=cols[j],lwd=6)}
	}

	n=dim(ALL2)[2]-1
	for (j in 1:11){
		selMale=which(sex=="Male" & colnames(ALL2)!=colnames(ALL2)[i])
		selFemale=which(sex=="Female" & colnames(ALL2)!=colnames(ALL2)[i])
		points(rep(j,length(selMale))+runif(length(selMale),-0.2,0.2),ALL2[j,selMale],pch=16,cex=0.4)
		points(rep(j,length(selFemale))+runif(length(selFemale),-0.2,0.2),ALL2[j,selFemale],pch=16,cex=0.4,col=3)
	}
	axis(1, at=1:11, labels=rownames(ALL2),las=2,cex.axis=0.7)

	text(beg,m*0.98,"Inferred sex:",cex=1,adj=0)
	text(beg,m*0.93,sex[i],cex=0.8,adj=0)
	legend(beg,m*0.8,c("No deletion","Deletion","Males","Females","Very low cov. (<10)","Low cov. (<20)"),col=c(4,2,1,3,"#ffcccc","#fde7b5"),pch=c(15,15,16,16,15,15),pt.cex=c(1,1,0.5,0.5,1,1),bg="white",y.intersp=1.5,cex=0.7)
	text(11.8,m*0.05,"* LCR is usually not covered\n by WES and can have\n high variability.",cex=0.6,adj=0)

	dev.off()

	pdf(file=paste(dir,"/06_copy-number/",rownames(ALL)[i],".copy-number_probas.pdf",sep=""),height=7,width=14)
	par(mfrow=c(1,2))

	plot(1:5,ylim=c(0,1),xlim=c(0.5,5.5), xaxt = "n",xlab="",ylab="Probability",cex=0,main=paste("Probabilities for CN of OPN1LW",sep=""),cex.main=1.5,cex.lab=1.2)
	axis(1, at=1:5, labels=c("0","1","2","3",">3"),cex.axis=1)
	text(x=3, y=-0.2, 'Full copy-number', xpd=NA,adj=0.5,cex=1.2)
	for (j in 1:4){
		rect(j-0.3,0,j+0.3,pLW[i,j],col=4)
	}
	rect(5-0.3,0,5+0.3,sum(pLW[i,5:8]),col=4)
	abline(h=0)

	plot(1:7,ylim=c(0,1),xlim=c(0.5,7.5), xaxt = "n",xlab="",ylab="Probability",cex=0,main=paste("Probabilities for CN of OPN1MW",sep=""),cex.main=1.5,cex.lab=1.2)
	axis(1, at=1:7, labels=c("0","1","2","3","4","5",">5"),cex.axis=1)
	text(x=4, y=-0.2, 'Full copy-number', xpd=NA,adj=0.5,cex=1.2)
	for (j in 1:6){
		rect(j-0.3,0,j+0.3,pMW[i,j],col=4)
	}
	rect(7-0.3,0,7+0.3,sum(pMW[i,7:10]),col=4)
	abline(h=0)

	dev.off()



}


if(dim(events)[1]>0){
	for (i in 1:dim(events)[1]){
		events[i,3]=gsub("\n"," ",events[i,3])
	}
}



write.table(events,file=paste(dir,"/0.events.tsv",sep=""),quote=F,row.names=F,sep="\t")
write.table(phenotype[,c(1,2,4,5,6,7)],file=paste(dir,"/0.phenotypes-all.tsv",sep=""),quote=F,row.names=F,sep="\t")
