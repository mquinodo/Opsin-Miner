
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
file9=args[9]
malesonly=args[10]
sexfile=args[11]

load(file6)

## reading files
# coverage from 0.coverage-ALL.tsv
data=read.table(file=file1,header=T,sep="\t",check.names = FALSE)
# protein haplotypes from 0.PROT.all.tsv
haplo=read.table(file=file2,skip=1)
# protein pathogenic haplotypes from 0.PROT.pathogenic.tsv
haploP=read.table(file=file3,skip=1)
haploD=haploP[,c(1,2,4,5,6)]
haploP=haploP[,c(1,3,4,5,6)]
# DNA haplotypes from 0.DNA.all.tsv
haploDNA=read.table(file=file4,skip=1)
# rare variants from 0.variants-rare.tsv
rareVariants=read.table(file=file5,header=T,row.names=NULL,sep="\t",colClasses = c("character"))
# coverage for Cys203Arg from 0.coverage.mut.ALL.tsv
muts=read.table(file=file8,skip=1)
# all variants from 0.variants-all.tsv
allVariants=read.table(file=file9,header=T,row.names=NULL,sep="\t",colClasses = c("character"))
# sexfile if provided
if(sexfile!="No"){
	sexinfo=read.table(file=sexfile,header=F,sep="\t")
}

# adding HGVS notation to the variants (from Variantvalidator)
allVariants=cbind(allVariants,allVariants[,1:3])
colnames(allVariants)[9:11]=c("HGVS-g.","HGVS-c.","HGVS-p.")
for (i in 1:dim(allVariants)[1]){
	sel=which(HGVS[,7]==allVariants[i,3] & HGVS[,8]==allVariants[i,4] & HGVS[,9]==allVariants[i,5])
	if(length(sel)==1){
		allVariants[i,9:11]=HGVS[sel,c(4,2,3)]
	} else {
		sel2=which(gnomAD2[,4]==allVariants[i,3] & gnomAD2[,5]==allVariants[i,4] & gnomAD2[,6]==allVariants[i,5])
		if(length(sel2)==1){
			allVariants[i,9:11]=c("NA",gnomAD2[sel2,c(9,8)])
		} else {
			allVariants[i,9:11]=c("NA","NA","NA")
		}
	}
}

# adding frequency from gnomAD and allelic ratio
allVariants=cbind(allVariants,allVariants[,1:7])
colnames(allVariants)[12:18]=c("AF-max-gnomAD","AF-OPN1LW-gnomAD","AF-OPN1MW-gnomAD","AF-OPN1MW2-gnomAD","Ratio","AF-gnomAD","AF-batch")
allVariants[,12]=0
allVariants[,13]=0
allVariants[,14]=0
allVariants[,15]=0
allVariants[,16]=0
allVariants[,17]=0
allVariants[,18]=0

if(dim(allVariants)[1]>0){
	for (j in 1:dim(allVariants)[1]){
		mut=as.numeric(allVariants[j,7])
		WT=as.numeric(allVariants[j,6])
		allVariants[j,16]=signif(mut/(mut+WT),digits=2)

		c=strsplit(allVariants[j,10],":")[[1]][2]
		FL=which(gnomAD[,7]==c & gnomAD[,1]=="OPN1LW")
		if(length(FL)==1){allVariants[j,13]=signif(gnomAD[FL,10],digits=2)}
		FL=which(gnomAD[,7]==c & gnomAD[,1]=="OPN1MW")
		if(length(FL)==1){allVariants[j,14]=signif(gnomAD[FL,10],digits=2)}
		FL=which(gnomAD[,7]==c & gnomAD[,1]=="OPN1MW2")
		if(length(FL)==1){allVariants[j,15]=signif(gnomAD[FL,10],digits=2)}
		allVariants[j,12]=max(allVariants[j,13:15])

		FL=which(gnomAD2[,4]==allVariants[j,3] & gnomAD2[,5]==allVariants[j,4] & gnomAD2[,6]==allVariants[j,5])
		if(length(FL)==1){allVariants[j,17]=signif(gnomAD2[FL,13],digits=2)}

		FL=which(allVariants[,3]==allVariants[j,3] & allVariants[,4]==allVariants[j,4] & allVariants[,5]==allVariants[j,5])
		allVariants[j,18]=signif(length(FL)/(dim(data)[2]-1),digits=2)
	}
}
allVariants=allVariants[,c(1:7,16,9:11,17,18)] # reordering columns

write.table(allVariants,file=paste(dir,"/03_variants/0.variants.annotated.tsv",sep=""),quote=F,row.names=F,sep="\t")
write.table(allVariants[which(allVariants[,12]<0.001 & allVariants[,13]<0.1),],file=paste(dir,"/03_variants/0.variants.annotated.rare.tsv",sep=""),quote=F,row.names=F,sep="\t")

# adding HGVS notation to the rare variants (from Variantvalidator)
rareVariants=cbind(rareVariants,rareVariants[,1:3])
colnames(rareVariants)[9:11]=c("HGVS-g.","HGVS-c.","HGVS-p.")
for (i in 1:dim(rareVariants)[1]){
	sel=which(HGVS[,7]==rareVariants[i,3] & HGVS[,8]==rareVariants[i,4] & HGVS[,9]==rareVariants[i,5])
	if(length(sel)==1){
		rareVariants[i,9:11]=HGVS[sel,c(4,2,3)]
	} else {
		sel2=which(gnomAD2[,4]==rareVariants[i,3] & gnomAD2[,5]==rareVariants[i,4] & gnomAD2[,6]==rareVariants[i,5])
		if(length(sel2)==1){
			rareVariants[i,9:11]=c("NA",gnomAD2[sel2,c(9,8)])
		} else {
			rareVariants[i,9:11]=c("NA","NA","NA")
		}
	}
}
if(length(which(rareVariants[,9]=="NA"))>0) {rareVariants=rareVariants[which(rareVariants[,9]!="NA"),]}

# normalize the coverage by control-chrX and remultiply by mean of control-chrX
data2=data
for (i in 2:dim(data2)[2]){
	data2[3:26,i]=as.numeric(data[3:26,i])/as.numeric(data[2,i])
}
data3=data2[,2:dim(data2)[2]]
rownames(data3)=data2[,1]
for (i in 1:dim(data3)[2]){
	data3[,i]=as.numeric(data3[,i])*mean(as.numeric(data[2,2:dim(data)[2]]))
}

# extracing coverage info for various targets (LCR and exons)
LW=1:dim(data3)[2]
MW=1:dim(data3)[2]
exons<-matrix(nrow=11,ncol=dim(data3)[2])
rownames(exons)=c("LCR*","Exon 3","Exon 6","Exon 1","Exon 2","Exon 4","Exon 5","Exon 1","Exon 2","Exon 4","Exon 5")
namesexons=c("LCR","ALL-exon3","ALL-exon6","OPN1LW-exon1","OPN1LW-exon2","OPN1LW-exon4","OPN1LW-exon5","OPN1MW-exon1","OPN1MW-exon2","OPN1MW-exon4","OPN1MW-exon5")
exons2=exons

# mean of normalized coverage per target without sample with less than 10% of average
me=data3[,1]
for (i in 1:dim(data3)[1]){
	me[i]=mean(as.numeric(data3[i,which(as.numeric(data3[i,])>0.1*mean(as.numeric(data3[i,])))]))
}

# normalize target coverage by mean and taking maximum for those exons covered by more than one target (all but LCR) -> exons
# taking maximum for those exons covered by more than one target (all but LCR) -> exons2 (no normalization)
data4=data[,-1]
for (i in 1:dim(data3)[2]){
	LW[i]=sum(data3[9:16,i])/sum(me[9:16]) # sum of coverage for exons specific to LW, normalized to average
	MW[i]=sum(data3[17:26,i])/sum(me[17:26]) # sum of coverage for exons specific to MW, normalized to average

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

# taking median values of raw values for LW and MW specific (medi1) and all targets (medi2)
medi1=c(median(apply(data4[9:16,],2,sum)),median(apply(data4[17:26,],2,sum)))
medi2=c(median(as.numeric(data4[3,])),median(apply(data4[5:6,],2,max)),median(apply(data4[7:8,],2,max)),
	median(apply(data4[9:10,],2,max)),median(apply(data4[11:12,],2,max)),median(apply(data4[13:14,],2,max)),median(apply(data4[15:16,],2,max)),
	median(apply(data4[21:22,],2,max)),median(apply(data4[23:24,],2,max)),median(apply(data4[25:26,],2,max)),median(apply(data4[17:18,],2,max)+apply(data4[19:20,],2,max)))

# change order to have "LCR","OPN1LW-exon1","OPN1LW-exon2","OPN1LW-exon4","OPN1LW-exon5","OPN1MW-exon1","OPN1MW-exon2","OPN1MW-exon4","OPN1MW-exon5","ALL-exon3","ALL-exon6"
exons=exons[c(1,4:11,2:3),]
exons2=exons2[c(1,4:11,2:3),]
medi2=medi2[c(1,4:11,2:3)]
namesexons=namesexons[c(1,4:11,2:3)]

## inferring sex
# ratio taken as control-chrX / control-autosome
rat=as.numeric(data[2,2:dim(data)[2]])/as.numeric(data[1,2:dim(data)[2]])

if(sexfile=="No"){
	# ratio for non large outliers
	rat2=rat[which(rat>mean(rat)-3*sd(rat) & rat<mean(rat)+3*sd(rat))]
	# taking 100 equally distant values between minimum and maximum ratio
	a=seq(min(rat2,na.rm=T)*1.1,max(rat2,na.rm=T)*0.9,(max(rat2,na.rm=T)*0.9-min(rat2,na.rm=T)*1.1)/100)
	# testing of assigning average to male and female samples among these 100 values
	b<-matrix(100000,ncol=101,nrow=101)
	for (i in 1:length(a)){
		for (j in 1:length(a)){
			if(i<j){
				a1=rat2[which(rat2>(a[i]+a[j])/2)] # presumed males
				a2=rat2[which(rat2<(a[i]+a[j])/2)] # presumed females
				b[i,j]=sum(abs(a1-a[j])^2)+sum(abs(a2-a[i])^2) # sum of deviations of male samples from tested male mean + same for females
			}
		}
	}
	# taking value that minimize the sum of deviations
	c=which(b == min(b,na.rm=T), arr.ind = TRUE)

	# now assigning sex to each sample using Z-score from the males / females distribution
	zm=rat
	zf=rat
	lim1=min(a[c[1]],a[c[2]])
	lim2=max(a[c[1]],a[c[2]])
	m=which(rat<(3*lim1+lim2)/4 & rat>(5*lim1-lim2)/4) # confident males 
	f=which(rat<(5*lim2-lim1)/4 & rat>(3*lim2+lim1)/4) # confident females
	# if one category empty, relax the criteria
	if(length(m)==0 | length(f)==0){
		m=which(rat<lim1+(lim2-lim1)/2 & rat>lim1-(lim2-lim1)/2) # confident males 
		f=which(rat<lim2+(lim2-lim1)/2 & rat>lim2-(lim2-lim1)/2) # confident females
	}
	for (i in 1:length(rat)){
		zm[i]=(rat[i]-mean(rat[m]))/sd(rat[m])
		zf[i]=(rat[i]-mean(rat[f]))/sd(rat[f])
		if(length(m)<2){
			zm[i]=(rat[i]-mean(rat[m]))/sd(rat[f])
		}
		if(length(f)<2){
			zf[i]=(rat[i]-mean(rat[f]))/sd(rat[m])
		}
	}
	z1=zm+zf
	sex=z1
	sex[which(z1<(-1.5))]="Male"
	sex[which(z1>(1.5))]="Female"
	sex[which(z1<=(1.5) & z1>=(-1.5))]="Undetermined"

	# if maleonly option is used, assign male to all samples
	if(malesonly=="Yes"){
		sex[1:length(sex)]="Male"
		sexout=cbind(colnames(data3),sex)
		colnames(sexout)=c("ID","Sex")
		write.table(sexout,file=paste(dir,"/04_coverage/0.inferred-sex.tsv",sep=""),quote=F,row.names=F,sep="\t")
	} else {
		sexout=cbind(colnames(data3),sex,z1)
		colnames(sexout)=c("ID","Inferred-sex","Z-score")
		write.table(sexout,file=paste(dir,"/04_coverage/0.inferred-sex.tsv",sep=""),quote=F,row.names=F,sep="\t")
		# printing figure with Z-score by sex
		pdf(file=paste(dir,"/04_coverage/","0.inferred-sex_Z-score.pdf",sep=""),height=10,width=20)
		plot(z1,col="white")
		points(which(z1<(-1.5)),z1[which(z1<(-1.5))],col=1)
		points(which(z1>(-1.5)),z1[which(z1>(-1.5))],col=3)
		points(which(z1<=(1.5) & z1>=(-1.5)),z1[which(z1<=(1.5) & z1>=(-1.5))],col=2)
		legend("topright",legend=c("Males","Females","Undetermined"),pch=1,col=c(1,3,2))
		dev.off()
		if(length(which(z1<(-1.5)))<3){print("Less than 3 males detected, please rerun with more males."); stop()}
		if(length(which(z1>1.5))<3){print("Less than 3 females detected, sex prediction is not accurate. Consider using --malesonly option if there are only males.")}
	}
} else {
	sex=2:dim(data)[2]
	for (i in 1:(dim(data)[2]-1)){
		if(length(which(sexinfo[,1]==colnames(data)[i+1]))==0){
			print(paste("Sex missing for sample ",colnames(data)[i+1],". Please add the information or run without --sexfile option.",sep=""))
			stop()
		}
		sex[i]=sexinfo[which(sexinfo[,1]==colnames(data)[i+1]),2]
	}
	if(length(which(sex=="Male"))<3){print("Less than 3 males in the input, please rerun with more males."); stop()}
}

# prepare out matrices
events<-matrix(ncol=9,nrow=0)
colnames(events)=c("#ID","Inferred-sex","Event","Pathogenicity","Details","Supporting-reads","WT-reads","Ratio","gnomAD-AF")
phenotype<-matrix(ncol=8,nrow=0)
colnames(phenotype)=c("#ID","Inferred-sex","Score","Phenotype","CN_OPN1LW","CN_OPN1MW","CN_hybrid","CN_hybrid2")

d=data[,-1]
ALL<-matrix(ncol=14,nrow=dim(d)[2])
rownames(ALL)=colnames(data)[2:dim(data)[2]]
colnames(ALL)=c("Control-autosome","Control-chrX","LCR","LCR-small","Exon3","Exon6","LWExon1","LWExon2","LWExon4","LWExon5","MWExon1","MWExon2","MWExon4","MWExon5")
# taking maximum for exons with two targets and normalize by control-autosome
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

# selecting male samples with at least 10 reads in all LW specific exons
LWmin=as.numeric(apply(ALL[,7:10],1,min)*d[1,])
sel=which(LWmin>=10 & sex=="Male")
# normalize all LW exons based on selected males and take minimum then renormalize by male samples
# this is the first estimate of the number of full LW genes
LWcor=cbind(ALL[,7]/mean(ALL[sel,7]),ALL[,8]/mean(ALL[sel,8]),ALL[,9]/mean(ALL[sel,9]),ALL[,10]/mean(ALL[sel,10]))
LWs=apply(LWcor,1,min)
LWs=LWs/mean(LWs[sel])

# same for MW
MWmin=as.numeric(apply(ALL[,11:14],1,min))
sel=which(MWmin>=mean(MWmin)/10 & sex=="Male")
MWcor=cbind(ALL[,11]/mean(ALL[sel,11]),ALL[,12]/mean(ALL[sel,12]),ALL[,13]/mean(ALL[sel,13]),ALL[,14]/mean(ALL[sel,14]))
MWs=apply(MWcor,1,min)
MWs=MWs/mean(MWs[sel])

# taking values between median*1.1 and median*3.1 and computing which value minimize devation from entire numbers (1,2,3,...,10)
f=seq(median(MWs[sel])*1.1,median(MWs[sel])*3.1,0.01)
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
# correct by value minimize deviation
MWs=MWs*f[which.min(r)]

# round estimates of number of copies and compute deviation
cLW=round(LWs)
cMW=round(MWs)
sdLW=sd(LWs-cLW)
sdMW=sd(MWs-cMW)

# compute probability for each copy number and each sample for LW
LWmax=8
pLW<-matrix(nrow=length(cLW),ncol=LWmax)
for (i in 1:length(cLW)){
    tot=1:LWmax
    # for 0 copies, SD would be 0 so we take all SD/4 instead
    tot[1]=dnorm(LWs[i],mean=0,sd=sdLW/4)
    for (j in 2:LWmax){
       tot[j]=dnorm(LWs[i],mean=j-1,sd=sdLW*j)
    }
    for (j in 1:LWmax){
       pLW[i,j]=tot[j]/sum(tot)
    }
}
# rounding with 2 digits
for (i in 1:dim(pLW)[2]){
	pLW[,i]=round(pLW[,i],digits=2)
}
# copy number with highest probability
cLW=apply(pLW,1,which.max)-1

# same for MW
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
for (i in 1:dim(pMW)[2]){
	pMW[,i]=round(pMW[,i],digits=2)
}

# putting in matrix and writing to text file
colnames(pLW)=c("0 copy","1 copy","2 copies","3 copies","4 copies","5 copies","6 copies","7 copies")
colnames(pMW)=c("0 copy","1 copy","2 copies","3 copies","4 copies","5 copies","6 copies","7 copies","8 copies","9 copies","10 copies","11 copies")
x=as.matrix(rownames(ALL))
colnames(x)="ID"
write.table(cbind(x,pLW),file=paste(dir,"/06_copy-number/0.proba_CN-LW.tsv",sep=""),quote=F,row.names=F,sep="\t")
write.table(cbind(x,pMW),file=paste(dir,"/06_copy-number/0.proba_CN-MW.tsv",sep=""),quote=F,row.names=F,sep="\t")

# normalization of each target/exon separately
ALL2=ALL

# LCR
i=3
LWmin=as.numeric(ALL[,i])
sel=which(LWmin>=mean(LWmin)/10 & sex=="Male")
ALL2[,i]=ALL[,i]/mean(ALL[sel,i])

# for exons 3 and 6
sel=which(sex=="Male")
ntot=mean(cLW[sel])+mean(cMW[sel]) # mean total number of copies per sample (LW and MW)
for(i in 5:6){
	LWmin=as.numeric(ALL[,i])
	# selecting males with more than 10% of mean coverage and less than 400% of mean coverage
	sel=which(LWmin>=mean(LWmin)/10 & sex=="Male" & LWmin<=mean(LWmin)*4)

	ALL2[,i]=ALL[,i]/mean(ALL[sel,i])
	# values from mean total number of copies -0.4 and +0.6
	f=seq(ntot-0.6,ntot+0.6,0.01)
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
	ALL2[,i]=ALL2[,i]*f[which.min(r)]
}

# same for exons 1,2,4,5 of LW
sel=which(sex=="Male")
ntot=mean(cLW[sel])
for(i in 7:10){
	LWmin=as.numeric(ALL[,i])
	sel=which(LWmin>=mean(LWmin)/10 & sex=="Male" & LWmin<=mean(LWmin)*4)
	ALL2[,i]=ALL[,i]/mean(ALL[sel,i])
	f=seq(ntot-0.6,ntot+0.6,0.01)
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
	ALL2[,i]=ALL2[,i]*f[which.min(r)]
}

# same for exons 1,2,4,5 of MW
sel=which(sex=="Male")
ntot=mean(cMW[sel])
for(i in 11:14){
	LWmin=as.numeric(ALL[,i])
	sel=which(LWmin>=mean(LWmin)/10 & sex=="Male" & LWmin<=mean(LWmin)*4)
	ALL2[,i]=ALL[,i]/mean(ALL[sel,i])
	f=seq(ntot-0.6,ntot+1,0.01)
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
	ALL2[,i]=ALL2[,i]*f[which.min(r)]
}

# reorder columns
ALL2=t(ALL2[,c(3,7,8,9,10,11,12,13,14,5,6)])
rownames(ALL2)=rownames(exons)

# looping through the samples to finish analysis and produce plots
for (i in 1:dim(data3)[2]){

	print(paste("Processing sample: ",rownames(ALL)[i],sep=""))

	# initialize scores for phenotype inference
	score1=0
	score2=0
	score3=0
	score4=0
	score5=0
	score6=0
	scoreA=0
	scoreB=0

	# thresholds for allelic frequency of events
	thr1=0.1
	thr2=0.4
	thr3=0.9

	# looking at haplotypes
	nreads=sum(haplo[which(haplo[,1]==colnames(data3)[i]),3])
	# pathogenic haplotype for the sample with more than 1 read
	x=haploP[which(haploP[,1]==colnames(data3)[i] & as.numeric(haploP[,3])>1),2:3]
	np=length(which(haploP[,1]==colnames(data3)[i] & as.numeric(haploP[,3])>1))
	if(length(np)>0){
		if(np==1){ # if 1 pathogenic haplotype
			# computing percentage of reads with pathogenic haplotype, updating score for phenotype inference and adding to events matrix
			perc=round(as.numeric(x[2])/(nreads),digits=2)
			pa="Low"
			if(as.numeric(x[2])>1 & as.numeric(perc)>thr3){
				pa="High"
				score1=score1+1
			}
			if(as.numeric(x[2])>1 & as.numeric(perc)<=thr3 & as.numeric(perc)>thr2){
				score2=score2+1
				scoreA=as.numeric(perc)
			}
			if(as.numeric(x[2])>1 & as.numeric(perc)<=thr2 & as.numeric(perc)>thr1){
				score3=score3+1
			}
			events=rbind(events,c(colnames(data3)[i],sex[i],"Pathogenic haplotype",pa,x[1,1],x[1,2],nreads-as.numeric(x[2]),perc,"NA"))
		}
		# same if more than 1 pathogenic haplotype is detected (very rare)
		if(np>1){
			x=x[sort(as.numeric(x[,2]), index.return=TRUE, decreasing=TRUE)$ix,]
			for (j in 1:dim(x)[1]){
				perc=round(as.numeric(x[j,2])/(nreads),digits=2)
				pa="Low"
				if(as.numeric(x[j,2])>1 & as.numeric(perc)>thr3){
					pa="High"
					score1=score1+1
				}
				if(as.numeric(x[j,2])>1 & as.numeric(perc)<=thr3 & as.numeric(perc)>thr2){
					score2=score2+1
				}
				if(as.numeric(x[j,2])>1 & as.numeric(perc)<=thr2 & as.numeric(perc)>thr1){
					score3=score3+1
				}
				events=rbind(events,c(colnames(data3)[i],sex[i],"Pathogenic haplotype",pa,x[j,1],x[j,2],nreads-as.numeric(x[j,2]),perc,"NA"))
			}
		}

	}

	# looking at rare variants: first annotate with gnomAD frequency, then update scores for phenotype and adding variants in events
	NS=rareVariants[which(rareVariants[,1]==colnames(data3)[i]),]
	if(dim(NS)[1]>0){
		for (j in 1:dim(NS)[1]){
			gnoMAX=0
			mut=as.numeric(NS[j,7])
			WT=as.numeric(NS[j,6])
			if(NS[j,11]!="NP_064445.2:p.?"){
				c=strsplit(NS[j,10],":")[[1]][2]
				gno=which(gnomAD[,7]==c)
				# if variant is only in one gene
				if(length(gno)==1){
					gno2=paste("gnomAD (",gnomAD[gno,1],"): ",signif(gnomAD[gno,10],digits=2),sep="")
					gnoMAX=max(gnoMAX,gnomAD[gno,10])
				}
				# if variant is in more than one gene
				if(length(gno)>1){
					gno2=paste("gnomAD (",gnomAD[gno[1],1],"): ",signif(gnomAD[gno[1],10],digits=2),sep="")
					gnoMAX=max(gnoMAX,gnomAD[gno[1],10])
					for (k in 2:length(gno)){
						gno2=paste(gno2,", gnomAD (",gnomAD[gno[k],1],"): ",signif(gnomAD[gno[k],10],digits=2),sep="")
						gnoMAX=max(gnoMAX,gnomAD[gno[k],10])
					}
				}
				if(length(gno)==0){gno2="not in gnomAD"} # if not in gnomAD
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
			if(mut>1 & mut/(WT+mut)>thr3 & gnoMAX<0.0001){
				pa="High"
			}
			# known pathogenic variants for BCM with high allelic frequency
			if(mut>1 & mut/(WT+mut)>thr3 & (grepl("p.\\(Cys203Arg\\)",NS[j,11])==T | grepl("p.\\(Asn94Lys\\)",NS[j,11])==T | grepl("p.\\(Arg330Gln\\)",NS[j,11])==T | grepl("p.\\(Gly338Glu\\)",NS[j,11])==T | grepl("p.\\(Trp177Arg\\)",NS[j,11])==T | grepl("Met1\\?",NS[j,11])==T  | grepl("Ter",NS[j,11])==T | grepl("\\+1G",NS[j,10])==T | grepl("\\+2T",NS[j,10])==T | grepl("-1G",NS[j,10])==T| grepl("-2A",NS[j,11])==T)){
				score1=score1+1
				pathoV=1
				scoreB=as.numeric(mut/(WT+mut))
			}
			if(mut>1 & grepl("p.\\(Trp177Arg\\)",NS[j,11])==T & mut/(WT+mut)>thr3){
				score6=score6+1
				scoreB=as.numeric(mut/(WT+mut))
			}
			if(mut>1 & mut/(WT+mut)>thr2 & mut/(WT+mut)<=thr3 & (grepl("p.\\(Cys203Arg\\)",NS[j,11])==T | grepl("p.\\(Asn94Lys\\)",NS[j,11])==T | grepl("p.\\(Arg330Gln\\)",NS[j,11])==T | grepl("p.\\(Gly338Glu\\)",NS[j,11])==T | grepl("p.\\(Trp177Arg\\)",NS[j,11])==T | grepl("Met1\\?",NS[j,11])==T  | grepl("Ter",NS[j,11])==T | grepl("\\+1G",NS[j,10])==T | grepl("\\+2T",NS[j,10])==T | grepl("-1G",NS[j,10])==T| grepl("-2A",NS[j,11])==T)){
				score2=score2+1
				pathoV=1
				scoreB=as.numeric(mut/(WT+mut))
			}
			if(mut>1 & mut/(WT+mut)>thr1 & mut/(WT+mut)<=thr2 & (grepl("p.\\(Cys203Arg\\)",NS[j,11])==T | grepl("p.\\(Asn94Lys\\)",NS[j,11])==T | grepl("p.\\(Arg330Gln\\)",NS[j,11])==T | grepl("p.\\(Gly338Glu\\)",NS[j,11])==T | grepl("p.\\(Trp177Arg\\)",NS[j,11])==T | grepl("Met1\\?",NS[j,11])==T  | grepl("Ter",NS[j,11])==T | grepl("\\+1G",NS[j,10])==T | grepl("\\+2T",NS[j,10])==T | grepl("-1G",NS[j,10])==T| grepl("-2A",NS[j,11])==T)){
				score3=score3+1
				pathoV=1
				scoreB=as.numeric(mut/(WT+mut))
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

	# looking at deletions if any coverage if less than 0.1 copy
	if(min(exons[,i])<=0.1){
		sco=0
		# LCR (if median is more than 20, it should be less than (1-thr3), or if median is more than 5, it should be 0 reads to be detected as deletion)
		if((medi2[1]>20 & exons[1,i]<(1-thr3)) | (medi2[1]>5 & exons[1,i]==0)){
			sco=max(sco,1)
		}
		score1=score1+sco
		# same exons 1,2,4,5 of LW
		sco=0
		for (k in 2:5){
			if((medi2[k]>20 & exons[k,i]<(1-thr3)) | (medi2[k]>5 & exons[k,i]==0)){
				sco=max(sco,1)
			}
		}
		score5=score5+sco
		# same exons 1,2,4,5 of MW
		sco=0
		for (k in 6:9){
			if((medi2[k]>20 & exons[k,i]<(1-thr3)) | (medi2[k]>5 & exons[k,i]==0)){
				sco=max(sco,1)
			}
		}
		score5=score5+sco
		# same for exons 3 and 6
		sco=0
		for (k in 10:11){
			if((medi2[k]>20 & exons[k,i]<(1-thr3)) | (medi2[k]>5 & exons[k,i]==0)){
				sco=max(sco,1)
			}
		}
		score1=score1+sco

		# adding deletions to events
		for (j in c(1:11)){
			if(exons[j,i]<0.1){
				pa="Low"
				if((medi2[j]>20 & exons[j,i]<(1-thr3)) | (medi2[j]>5 & exons[j,i]==0)){
					pa="High"
				}
				events=rbind(events,c(colnames(data3)[i],sex[i],"Hemi/homozygous deletion",pa,namesexons[j],"NA",exons2[j,i],round(exons[j,i],digits=2),"NA"))
			}
		}
	}

	# starting PDF and first plot
	pdf(file=paste(dir,"/07_plots/",rownames(ALL)[i],"_OPN1LW-MW_analysis.pdf",sep=""),height=7,width=14)
	par(mfrow=c(1,2),mar=c(2, 4, 4, 2) + 0.1)
	plot(1,xlim=c(0,10),ylim=c(0,9.5),xaxt ="n",yaxt ="n",xlab="",ylab="",cex=0,main=paste("Summary for ",colnames(ALL2)[i],sep=""),cex.main=1.5)
	par(xpd=NA)
	par(xpd=FALSE)
	text(0,9.3,"Graphical representation and inferred copy-numbers:",cex=1,adj=0)
	m=8.2
	lines(c(0,10),c(m,m))
	# plotting graphical representation of LCR
	if((medi2[1]>20 & exons[1,i]<0.1) | (medi2[1]>5 & exons[1,i]==0)){
		rect(0.2,m-0.2,1,m+0.2,col="white",border="NA")
		lines(c(0.25,0.2,0.2,0.25),c(m-0.2,m-0.2,m+0.2,m+0.2))
		lines(c(0.95,1,1,0.95),c(m-0.2,m-0.2,m+0.2,m+0.2))
		text(0.6,m+0.5,"LCR deleted",col="red",cex=0.8)
	} else {
		rect(0.2,m-0.2,1,m+0.2,col="#7c7c7d")
		text(0.6,m+0.5,"LCR",cex=0.8)
	}

	# number of copies of LW (+1!, so 2 means one copy)
	nLW=which.max(pLW[i,])
	# correction for LW if hybrids (if one is <0.8*mean and has 3 copies detected or if one exon is less than 0.1 and more than 1 copies are detected)
	if(min(ALL2[2:5,i])<0.8*mean(ALL2[2:5,i]) & nLW==4){nLW=3}
	if(min(ALL2[2:5,i])<0.1 & nLW>1){nLW=1}	
	# number of copies of MW (+1!)
	nMW=which.max(pMW[i,])
	# if low coverage of exon 3 or 6 then copies of LW and MW should be 0 (+1)
	if(ALL[i,5]<mean(ALL[,5])/10 | ALL[i,6]<mean(ALL[,6])/10){
		nLW=1
		nMW=1
	}

	# graphical parameters for plotting
	xl=7.3
	hybrid=0
	hybrid2=0
	hybrid3=0
	CH=0
	CH2=0
	CH3=0
	e1=0.292
	g1=0.15
	l1=0.3
	l2=0.2
	y1=m-l1/2
	y2=y1+l1
	x1=2
	colo3=c()
	colo4=c()
	fact=1.2
	poshyb=0

	# checking if there could be a full hybrid gene
	# checking that exon3 / exon6 / exon1 / exon2 / ... are at least covered in one gene
	if(ALL[i,5]>mean(ALL[,5])/10 & ALL[i,6]>mean(ALL[,6])/10 & (ALL[i,7]>mean(ALL[,7])/10 | ALL[i,11]>mean(ALL[,11])/10) & (ALL[i,8]>mean(ALL[,8])/10 | ALL[i,12]>mean(ALL[,12])/10) & (ALL[i,9]>mean(ALL[,9])/10 | ALL[i,13]>mean(ALL[,13])/10) & (ALL[i,10]>mean(ALL[,10])/10 | ALL[i,14]>mean(ALL[,14])/10) ){
		# checking that there is at least one exon covered in LW and in MW specific exons
		if((ALL[i,7]>mean(ALL[,7])/10 | ALL[i,8]>mean(ALL[,8])/10 | ALL[i,9]>mean(ALL[,9])/10 | ALL[i,10]>mean(ALL[,10])/10) & (ALL[i,11]>mean(ALL[,11])/10 | ALL[i,12]>mean(ALL[,12])/10 | ALL[i,13]>mean(ALL[,13])/10 | ALL[i,14]>mean(ALL[,14])/10)){
			poshyb=1
		}
	}

	# if there is no full LW gene, at least one full MW and the precence of hybrid is possible
	if(nLW==1 & nMW>1 & poshyb==1){
		# initialize colors of exons as red (LW) and grey for 3 and 6
		colo=c("#ff6f69","#ff6f69","grey","#ff6f69","#ff6f69","grey")
		nMWtest=min(which(pMW[i,]>0.2))
		if(nMWtest>5){nMWtest=nMWtest-1} # correction for high number of MW
		if(ALL[i,7]<mean(ALL[,7])/10 & ALL2[6,i]>fact*(nMWtest-1)){colo[1]="#88d8b0"} # if exon 1 of LW is deleted and exon 1 of MW has coverage higher than (fact*expected copy number for MW)
		if(ALL[i,8]<mean(ALL[,8])/10 & ALL2[7,i]>fact*(nMWtest-1)){colo[2]="#88d8b0"} # same for exon 2
		if(ALL[i,9]<mean(ALL[,9])/10 & ALL2[8,i]>fact*(nMWtest-1)){colo[4]="#88d8b0"} # same for exon 4
		if(ALL[i,10]<mean(ALL[,10])/10 & ALL2[9,i]>fact*(nMWtest-1)){colo[5]="#88d8b0"} # same for exon 5
		# plotting graphical representation of the hybrid gene
		if(colo[1]!="#ff6f69" | colo[2]!="#ff6f69" | colo[4]!="#ff6f69" | colo[5]!="#ff6f69"){
			hybrid=1
			for(d in 1:6){
				rect(x1+e1*(d-1)+g1*(d-1),y1,x1+e1*d+g1*(d-1),y2,col=colo[d])
				text(x1+e1*(d-1)+g1*(d-1)+e1/2,m,d,cex=0.8)
			}
			# computing copy number of the hybrid gene
			ok=c()
			if(ALL[i,7]>mean(ALL[,7])/10){ok=c(ok,2)}
			if(ALL[i,8]>mean(ALL[,8])/10){ok=c(ok,3)}
			if(ALL[i,9]>mean(ALL[,9])/10){ok=c(ok,4)}
			if(ALL[i,10]>mean(ALL[,10])/10){ok=c(ok,5)}
			CH=round(mean(ALL2[ok,i]))
			# adding hybrid in events
			colo2=colo
			colo2[which(colo2=="#ff6f69")]="L"
			colo2[which(colo2=="#88d8b0")]="M"
			colo2[which(colo2=="grey")]="x"
			colo3=paste(colo2,collapse="")
			events=rbind(events,c(colnames(data3)[i],sex[i],"Hybrid gene","High",colo3,"NA","NA","NA","NA"))
		}
	}
	# same if there is no full MW gene, at least one full LW and the precence of hybrid is possible
	if(nLW>1 & nMW==1 & poshyb==1){
		colo=c("#88d8b0","#88d8b0","grey","#88d8b0","#88d8b0","grey")
		nLWtest=min(which(pLW[i,]>0.2))
		if(ALL[i,11]<mean(ALL[,11])/10 & ALL2[2,i]>fact*(nLWtest-1)){colo[1]="#ff6f69"}
		if(ALL[i,12]<mean(ALL[,12])/10 & ALL2[3,i]>fact*(nLWtest-1)){colo[2]="#ff6f69"}
		if(ALL[i,13]<mean(ALL[,13])/10 & ALL2[4,i]>fact*(nLWtest-1)){colo[4]="#ff6f69"}
		if(ALL[i,14]<mean(ALL[,14])/10 & ALL2[5,i]>fact*(nLWtest-1)){colo[5]="#ff6f69"}
		if(colo[1]!="#88d8b0" | colo[2]!="#88d8b0" | colo[4]!="#88d8b0" | colo[5]!="#88d8b0"){
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
			colo2[which(colo2=="#ff6f69")]="L"
			colo2[which(colo2=="#88d8b0")]="M"
			colo2[which(colo2=="grey")]="x"
			colo3=paste(colo2,collapse="")
			events=rbind(events,c(colnames(data3)[i],sex[i],"Hybrid gene","High",colo3,"NA","NA","NA","NA"))
		}
	}
	# same if hybrid is possible but there are no full LW nor MW
	if(nLW==1 & nMW==1 & poshyb==1){
		colo=c("#ff6f69","#ff6f69","grey","#ff6f69","#ff6f69","grey")
		if(ALL[i,7]<mean(ALL[,7])/10){colo[1]="#88d8b0"}
		if(ALL[i,8]<mean(ALL[,8])/10){colo[2]="#88d8b0"}
		if(ALL[i,9]<mean(ALL[,9])/10){colo[4]="#88d8b0"}
		if(ALL[i,10]<mean(ALL[,10])/10){colo[5]="#88d8b0"}
		hybrid=1
		CH=1
		for(d in 1:6){
			rect(x1+e1*(d-1)+g1*(d-1),y1,x1+e1*d+g1*(d-1),y2,col=colo[d])
			text(x1+e1*(d-1)+g1*(d-1)+e1/2,m,d,cex=0.8)
		}
		colo2=colo
		colo2[which(colo2=="#ff6f69")]="L"
		colo2[which(colo2=="#88d8b0")]="M"
		colo2[which(colo2=="grey")]="x"
		colo3=paste(colo2,collapse="")
		events=rbind(events,c(colnames(data3)[i],sex[i],"Hybrid gene","High",colo3,"NA","NA","NA","NA"))
	}

	# computing if additional hybrids could exist
	limi=0.25
	# taking copy number of exons1-2-4-5 LW, exons1-2-4-5 MW, exons 3-6 and then remove copies from full LW and full MW genes
	initial=ALL2[2:11,i]
	initial[c(1:4,9:10)]=initial[c(1:4,9:10)]-(nLW-1)
	initial[c(5:8,9:10)]=initial[c(5:8,9:10)]-(nMW-1)
	if(poshyb==1){
		# if there is already an hybrid detected, we also substract relevant exons
		if(hybrid==1){
			if(initial[1]>limi & initial[2]>limi & initial[3]>limi & initial[8]>limi & initial[9]>limi & initial[10]>limi){
				initial[c(1,2,3,8,9,10)]=initial[c(1,2,3,8,9,10)]-1
			} else if (initial[1]>limi & initial[2]>limi & initial[7]>limi & initial[4]>limi & initial[9]>limi & initial[10]>limi) {
				initial[c(1,2,7,4,9,10)]=initial[c(1,2,7,4,9,10)]-1
			} else if (initial[1]>limi & initial[6]>limi & initial[3]>limi & initial[4]>limi & initial[9]>limi & initial[10]>limi) {
				initial[c(1,6,3,4,9,10)]=initial[c(1,6,3,4,9,10)]-1
			} else if (initial[5]>limi & initial[2]>limi & initial[3]>limi & initial[4]>limi & initial[9]>limi & initial[10]>limi) {
				initial[c(5,2,3,4,9,10)]=initial[c(5,2,3,4,9,10)]-1
			} else if (initial[1]>limi & initial[6]>limi & initial[7]>limi & initial[8]>limi & initial[9]>limi & initial[10]>limi) {
				initial[c(1,6,7,8,9,10)]=initial[c(1,6,7,8,9,10)]-1
			} else if (initial[5]>limi & initial[2]>limi & initial[7]>limi & initial[8]>limi & initial[9]>limi & initial[10]>limi) {
				initial[c(5,2,7,8,9,10)]=initial[c(5,2,7,8,9,10)]-1
			} else if (initial[5]>limi & initial[6]>limi & initial[3]>limi & initial[8]>limi & initial[9]>limi & initial[10]>limi) {
				initial[c(5,6,3,8,9,10)]=initial[c(5,6,3,8,9,10)]-1
			} else if (initial[5]>limi & initial[6]>limi & initial[7]>limi & initial[4]>limi & initial[9]>limi & initial[10]>limi) {
				initial[c(5,6,7,4,9,10)]=initial[c(5,6,7,4,9,10)]-1
			} else if (initial[1]>limi & initial[2]>limi & initial[7]>limi & initial[8]>limi & initial[9]>limi & initial[10]>limi) {
				initial[c(1,2,7,8,9,10)]=initial[c(1,2,7,8,9,10)]-1
			} else if (initial[5]>limi & initial[6]>limi & initial[3]>limi & initial[4]>limi & initial[9]>limi & initial[10]>limi) {
				initial[c(5,6,3,4,9,10)]=initial[c(5,6,3,4,9,10)]-1
			} else if (initial[1]>limi & initial[6]>limi & initial[7]>limi & initial[4]>limi & initial[9]>limi & initial[10]>limi) {
				initial[c(1,6,7,4,9,10)]=initial[c(1,6,7,4,9,10)]-1
			} else if (initial[5]>limi & initial[2]>limi & initial[3]>limi & initial[8]>limi & initial[9]>limi & initial[10]>limi) {
				initial[c(5,2,3,8,9,10)]=initial[c(5,2,3,8,9,10)]-1
			} else if (initial[1]>limi & initial[6]>limi & initial[3]>limi & initial[8]>limi & initial[9]>limi & initial[10]>limi) {
				initial[c(1,6,3,8,9,10)]=initial[c(1,6,3,8,9,10)]-1
			} else if (initial[5]>limi & initial[2]>limi & initial[7]>limi & initial[4]>limi & initial[9]>limi & initial[10]>limi) {
				initial[c(5,2,7,4,9,10)]=initial[c(5,2,7,4,9,10)]-1
			}
		}
		# if there is still signal from exons 1 to 6 (from any gene), we have an additional hybrid (hybrid2) and we add it in events
		if(max(initial[1],initial[2],initial[3],initial[4])>limi & max(initial[5],initial[6],initial[7],initial[8])>limi & max(initial[1],initial[5])>limi & max(initial[2],initial[6])>limi & max(initial[3],initial[7])>limi & max(initial[4],initial[8])>limi & initial[9]>limi & initial[10]>limi){
			hybrid2=1
			CH2=1
			colo2=c("L","L","x","L","L","x")
			if(initial[5]>initial[1]){colo2[1]="M"; initial[5]=initial[5]-1} else {initial[1]=initial[1]-1}
			if(initial[6]>initial[2]){colo2[2]="M"; initial[6]=initial[6]-1} else {initial[2]=initial[2]-1}
			if(initial[7]>initial[3]){colo2[4]="M"; initial[7]=initial[7]-1} else {initial[3]=initial[3]-1}
			if(initial[8]>initial[4]){colo2[5]="M"; initial[8]=initial[8]-1} else {initial[4]=initial[4]-1}
			initial[9]=initial[9]-1
			initial[10]=initial[10]-1
			colo3=paste(colo2,collapse="")

			# correction if not hybrid detected
			if(colo2[1]==colo2[2] & colo2[1]==colo2[4] & colo2[1]==colo2[5]){
				hybrid2=0
				CH2=0
			} else {
				events=rbind(events,c(colnames(data3)[i],sex[i],"Possible additional hybrid gene","Low",colo3,"NA","NA","NA","NA"))
			}
		}
		if(max(initial[1],initial[2],initial[3],initial[4])>limi & max(initial[5],initial[6],initial[7],initial[8])>limi & max(initial[1],initial[5])>limi & max(initial[2],initial[6])>limi & max(initial[3],initial[7])>limi & max(initial[4],initial[8])>limi & initial[9]>limi & initial[10]>limi){
			hybrid3=1
			CH3=1
			colo2=c("L","L","x","L","L","x")
			if(initial[5]>initial[1]){colo2[1]="M"; initial[5]=initial[5]-1} else {initial[1]=initial[1]-1}
			if(initial[6]>initial[2]){colo2[2]="M"; initial[6]=initial[6]-1} else {initial[2]=initial[2]-1}
			if(initial[7]>initial[3]){colo2[4]="M"; initial[7]=initial[7]-1} else {initial[3]=initial[3]-1}
			if(initial[8]>initial[4]){colo2[5]="M"; initial[8]=initial[8]-1} else {initial[4]=initial[4]-1}
			initial[9]=initial[9]-1
			initial[10]=initial[10]-1
			colo4=paste(colo2,collapse="")
			
			# correction if not hybrid detected
			if(colo2[1]==colo2[2] & colo2[1]==colo2[4] & colo2[1]==colo2[5]){
				hybrid3=0
				CH3=0
			} else {
				events=rbind(events,c(colnames(data3)[i],sex[i],"Possible additional hybrid gene","Low",colo4,"NA","NA","NA","NA"))
			}
		} 
	}

	# correction for hybrid
	corr=c(0,0,0,0)
	if(abs(initial[4]-initial[8])>1 & abs(initial[4])>0.5 & abs(initial[8])>0.5 & initial[4]>initial[8] & hybrid==0 & length(which(initial==0))==0 & nMW>2){corr[4]=1}
	if(abs(initial[4]-initial[8])>1 & abs(initial[4])>0.5 & abs(initial[8])>0.5 & initial[4]<initial[8] & hybrid==0 & length(which(initial==0))==0 & nLW>2){corr[4]=2}
	if(abs(initial[3]-initial[7])>1 & abs(initial[3])>0.5 & abs(initial[7])>0.5 & initial[3]>initial[7] & hybrid==0 & length(which(initial==0))==0 & nMW>2){corr[3]=1}
	if(abs(initial[3]-initial[7])>1 & abs(initial[3])>0.5 & abs(initial[7])>0.5 & initial[3]<initial[7] & hybrid==0 & length(which(initial==0))==0 & nLW>2){corr[3]=2}
	if(abs(initial[2]-initial[6])>1 & abs(initial[2])>0.5 & abs(initial[6])>0.5 & initial[2]>initial[6] & hybrid==0 & length(which(initial==0))==0 & nMW>2){corr[2]=1}
	if(abs(initial[2]-initial[6])>1 & abs(initial[2])>0.5 & abs(initial[6])>0.5 & initial[2]<initial[6] & hybrid==0 & length(which(initial==0))==0 & nLW>2){corr[2]=2}
	if(abs(initial[1]-initial[5])>1 & abs(initial[1])>0.5 & abs(initial[5])>0.5 & initial[1]>initial[5] & hybrid==0 & length(which(initial==0))==0 & nMW>2){corr[1]=1}
	if(abs(initial[1]-initial[5])>1 & abs(initial[1])>0.5 & abs(initial[5])>0.5 & initial[1]<initial[5] & hybrid==0 & length(which(initial==0))==0 & nLW>2){corr[1]=2}
	if(length(which(corr>0)==1)){
		if(corr[1]==2){nLW=nLW-1; initial[1]=initial[1]+1; initial[5]=initial[5]-1; hybrid2=1; CH2=1; colo3="MLxLLx"; events=rbind(events,c(colnames(data3)[i],sex[i],"Possible additional hybrid gene","Low",colo3,"NA","NA","NA","NA"))}
		if(corr[1]==1){nMW=nMW-1; initial[1]=initial[1]-1; initial[5]=initial[5]+1; hybrid2=1; CH2=1; colo3="LMxMMx"; events=rbind(events,c(colnames(data3)[i],sex[i],"Possible additional hybrid gene","Low",colo3,"NA","NA","NA","NA"))}
		if(corr[2]==2){nLW=nLW-1; initial[2]=initial[2]+1; initial[6]=initial[6]-1; hybrid2=1; CH2=1; colo3="LMxLLx"; events=rbind(events,c(colnames(data3)[i],sex[i],"Possible additional hybrid gene","Low",colo3,"NA","NA","NA","NA"))}
		if(corr[2]==1){nMW=nMW-1; initial[2]=initial[2]-1; initial[6]=initial[6]+1; hybrid2=1; CH2=1; colo3="MLxMMx"; events=rbind(events,c(colnames(data3)[i],sex[i],"Possible additional hybrid gene","Low",colo3,"NA","NA","NA","NA"))}
		if(corr[3]==2){nLW=nLW-1; initial[3]=initial[3]+1; initial[7]=initial[7]-1; hybrid2=1; CH2=1; colo3="LLxMLx"; events=rbind(events,c(colnames(data3)[i],sex[i],"Possible additional hybrid gene","Low",colo3,"NA","NA","NA","NA"))}
		if(corr[3]==1){nMW=nMW-1; initial[3]=initial[3]-1; initial[7]=initial[7]+1; hybrid2=1; CH2=1; colo3="MMxLMx"; events=rbind(events,c(colnames(data3)[i],sex[i],"Possible additional hybrid gene","Low",colo3,"NA","NA","NA","NA"))}
		if(corr[4]==2){nLW=nLW-1; initial[4]=initial[4]+1; initial[8]=initial[8]-1; hybrid2=1; CH2=1; colo3="LLxLMx"; events=rbind(events,c(colnames(data3)[i],sex[i],"Possible additional hybrid gene","Low",colo3,"NA","NA","NA","NA"))}
		if(corr[4]==1){nMW=nMW-1; initial[4]=initial[4]-1; initial[8]=initial[8]+1; hybrid2=1; CH2=1; colo3="MMxMLx"; events=rbind(events,c(colnames(data3)[i],sex[i],"Possible additional hybrid gene","Low",colo3,"NA","NA","NA","NA"))}
	}

	x1=2
	# plotting LW gene as a red arrow
	if(nLW>1){
		x2=4
		x3=4.5
		y1=m-l1/2
		y2=y1+l1
		y3=y1-l2
		y4=y2+l2
		y5=(y1+y2)/2
		polygon(c(x1,x1,x2,x2,x3,x2,x2),c(y1,y2,y2,y4,y5,y3,y1),col="#ff6f69",border="NA")
		lines(c(x1,x1,x2,x2,x3,x2,x2,x1),c(y1,y2,y2,y4,y5,y3,y1,y1))
		text(x1+(x3-x1)/2,m+0.5,"OPN1LW",cex=0.8)
		if(nLW<=2) {
			text(x1+(x3-x1)/2,m-0.5,paste(nLW-1," copy",sep=""),cex=0.8)
		} else {
			text(x1+(x3-x1)/2,m-0.5,paste(nLW-1," copies",sep=""),cex=0.8)
		}
	}
	# plotting if no copy of LW and no hybrid
	if(nLW==1 & hybrid==0){
		colo=c("#ff6f69","#ff6f69","grey","#ff6f69","#ff6f69","grey")
		nLWtest=min(which(pLW[i,]>0.2))
		nMWtest=min(which(pMW[i,]>0.2))
		# checking which exons are deleted and which not
		good=c(1:6)
		bad=c()
		if(ALL[i,7]<mean(ALL[,7])/10){bad=c(bad,1); good=good[-which(good==1)]}
		if(ALL[i,8]<mean(ALL[,8])/10){bad=c(bad,2); good=good[-which(good==2)]}
		if(ALL[i,9]<mean(ALL[,9])/10){bad=c(bad,4); good=good[-which(good==4)]}
		if(ALL[i,10]<mean(ALL[,10])/10){bad=c(bad,5); good=good[-which(good==5)]}
		if(ALL[i,10]<mean(ALL[,10])/10 & ALL2[11,i]-1<fact*(nLWtest+nMWtest-2)){bad=c(bad,6); good=good[-which(good==6)];}
		if(ALL[i,6]<mean(ALL[,6])/10){bad=c(bad,6); good=good[which(good!=6)];}
		if(ALL[i,8]<mean(ALL[,8])/10 | ALL[i,9]<mean(ALL[,9])/10){bad=c(bad,3); good=good[-which(good==3)]}
		# plotting
		for(d in bad){
			rect(x1+e1*(d-1)+g1*(d-1),y1,x1+e1*d+g1*(d-1),y2,col="white",border="NA")
			lines(c(x1+e1*(d-1)+g1*(d-1)+0.05,x1+e1*(d-1)+g1*(d-1),x1+e1*(d-1)+g1*(d-1),x1+e1*(d-1)+g1*(d-1)+0.05),c(m-0.2,m-0.2,m+0.2,m+0.2))
			lines(c(x1+e1*d+g1*(d-1)-0.05,x1+e1*d+g1*(d-1),x1+e1*d+g1*(d-1),x1+e1*d+g1*(d-1)-0.05),c(m-0.2,m-0.2,m+0.2,m+0.2))
		}
		for(d in good){
			rect(x1+e1*(d-1)+g1*(d-1),y1,x1+e1*d+g1*(d-1),y2,col=colo[d])
			text(x1+e1*(d-1)+g1*(d-1)+e1/2,m,d,cex=0.8)
		}
		# plotting in case of full LW deletion
		if(ALL[i,7]<mean(ALL[,7])/10 & ALL[i,8]<mean(ALL[,8])/10 & ALL[i,9]<mean(ALL[,9])/10 & ALL[i,10]<mean(ALL[,10])/10){
			rect(2,m-0.22,4.5,m+0.22,col="white",border="NA")
			lines(c(4.45,4.5,4.5,4.45),c(m-0.2,m-0.2,m+0.2,m+0.2))
			lines(c(2.05,2,2,2.05),c(m-0.2,m-0.2,m+0.2,m+0.2))
			text(3.25,m+0.5,"OPN1LW deleted",col="red",cex=0.8)
		}
		good=good[which(good!=3 & good!=6)]
		if(length(good)>0){
			text(3.25,m+0.5,"Partial OPN1LW",col="red",cex=0.8)
		}
	}
	# if no LW but hybrid is there (already plotted before)
	if(nLW==1 & hybrid==1){
		x1=2
		x2=4
		x3=4.5
		text(x1+(x3-x1)/2,m+0.5,"Hybrid gene",cex=0.8)
		if(CH>1){
			text(x1+(x3-x1)/2,m-0.5,paste(CH," copies",sep=""),cex=0.8)
		}
		if(CH==1){
			text(x1+(x3-x1)/2,m-0.5,paste(CH," copy",sep=""),cex=0.8)
		}
	}
	# if no MW and no hybrid, plot incomplete gene or deletion
	if(nMW==1 & hybrid==0){
		x1=6
		colo=c("#88d8b0","#88d8b0","grey","#88d8b0","#88d8b0","grey")
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
		if(ALL[i,11]<mean(ALL[,11])/10 & ALL[i,12]<mean(ALL[,12])/10 & ALL[i,13]<mean(ALL[,13])/10 & ALL[i,14]<mean(ALL[,14])/10){
			rect(6,m-0.22,8.5,m+0.22,col="white",border="NA")
			lines(c(8.45,8.5,8.5,8.45),c(m-0.2,m-0.2,m+0.2,m+0.2))
			lines(c(6.05,6,6,6.05),c(m-0.2,m-0.2,m+0.2,m+0.2))
			text(7.25,m+0.5,"OPN1MW deleted",col="red",cex=0.8)
		}
		good=good[which(good!=3 & good!=6)]
		if(length(good)>0){
			text(7.25,m+0.5,"Partial OPN1MW",col="red",cex=0.8)
		}
	}
	# if hybrid and LW (already plotted)
	if(nMW==1 & hybrid==1 & nLW!=1){
		x1=6
		x2=8
		x3=8.5
		text(x1+(x3-x1)/2,m+0.5,"Hybrid gene",cex=0.8)
		if(CH>1){
			text(x1+(x3-x1)/2,m-0.5,paste(CH," copies",sep=""),cex=0.8)
		}
		if(CH==1){
			text(x1+(x3-x1)/2,m-0.5,paste(CH," copy",sep=""),cex=0.8)
		}
	}
	# if 1 or more MW, plot as green arrow
	if(nMW>1){
		x1=6
		x2=8
		x3=8.5
		y1=m-l1/2
		y2=y1+l1
		y3=y1-l2
		y4=y2+l2
		y5=(y1+y2)/2
		polygon(c(x1,x1,x2,x2,x3,x2,x2),c(y1,y2,y2,y4,y5,y3,y1),col="#88d8b0",border="NA")
		lines(c(x1,x1,x2,x2,x3,x2,x2,x1),c(y1,y2,y2,y4,y5,y3,y1,y1))
		text(x1+(x3-x1)/2,m+0.5,"OPN1MW",cex=0.8)
		if(nMW<=2) {
			text(x1+(x3-x1)/2,m-0.5,paste(nMW-1," copy",sep=""),cex=0.8)
		} else {
			text(x1+(x3-x1)/2,m-0.5,paste(nMW-1," copies",sep=""),cex=0.8)
		}
	}

	# write if possible other hybrid genes have been found
	xl=7
	xlt=7
	if(hybrid2==1 & hybrid==0 & hybrid3==0){
		text(5,7.3,paste("Possible additional hybrid gene (",colo3,")",sep=""),cex=0.9,adj=0.5,col="red")
		xl=xl-0.2
		xlt=xlt-0.2
	}
	if(hybrid2==1 & hybrid==1 & hybrid3==0){
		text(5,7.3,paste("Possible other hybrid gene (",colo3,")",sep=""),cex=0.9,adj=0.5,col="red")
		xl=xl-0.2
		xlt=xlt-0.2
	}
	if(hybrid2==1 & hybrid==0 & hybrid3==1){
		text(5,7.3,paste("Possible additional hybrid genes (",colo3," and ",colo4,")",sep=""),cex=0.9,adj=0.5,col="red")
		xl=xl-0.2
		xlt=xlt-0.2
	}
	if(hybrid2==1 & hybrid==1 & hybrid3==1){
		text(5,7.3,paste("Possible other hybrid genes (",colo3," and ",colo4,")",sep=""),cex=0.9,adj=0.5,col="red")
		xl=xl-0.2
		xlt=xlt-0.2
	}

	# writing protein haplotypes on plot
	nreads=sum(haplo[which(haplo[,1]==colnames(data3)[i]),3])
	text(0,xl,"Exon 3 protein haplotypes:",adj=0,cex=1)
	x=haplo[which(haplo[,1]==colnames(data3)[i]),2:3]
	x=x[sort(as.numeric(x[,2]), index.return=TRUE, decreasing=TRUE)$ix,]

	pathoP=unique(haploP[,2])

	if(dim(x)[1]>0){
		for (j in 1:dim(x)[1]){
			perc=round(as.numeric(x[j,2])/(nreads),digits=2)
			if(is.element(x[j,1],pathoP)==T){
				if(perc>0.1 | x[j,2]>4){
					text(0,xl-0.4,paste(x[j,1]," (",perc*100,"%, n=",x[j,2],"/",nreads,", pathogenic)",sep=""),adj=0,cex=0.8,col="red")
					xl=xl-0.4
				}
			} else {
				if(perc>0.1 | x[j,2]>4){
					text(0,xl-0.4,paste(x[j,1]," (",perc*100,"%, n=",x[j,2],"/",nreads,")",sep=""),adj=0,cex=0.8)
					xl=xl-0.4
				}
			}		
		}
	} else {
		text(0,xl-0.4,"No haplotype detected",adj=0,cex=0.8)
		xl=xl-0.4
	}
	xl=xl+0.15

	# writing DNA haplotypes on plot
	nreads=sum(haploDNA[which(haploDNA[,1]==colnames(data3)[i]),3])
	text(5,xlt,"Exon 3 DNA haplotypes:",adj=0,cex=1)
	x=haploDNA[which(haploDNA[,1]==colnames(data3)[i]),2:3]
	x=x[sort(as.numeric(x[,2]), index.return=TRUE, decreasing=TRUE)$ix,]

	pathoD=unique(haploD[,2])

	if(dim(x)[1]>0){
		for (j in 1:dim(x)[1]){
			perc=round(as.numeric(x[j,2])/(nreads),digits=2)
			if(is.element(x[j,1],pathoD)==T){
				if(perc>0.1 | x[j,2]>4){
					text(5,xlt-0.4,paste(x[j,1]," (",perc*100,"%, n=",x[j,2],"/",nreads,", pathogenic)",sep=""),adj=0,cex=0.7,col="red")
					xlt=xlt-0.3
				}
			} else {	
				if(perc>0.1 | x[j,2]>4){
					text(5,xlt-0.4,paste(x[j,1]," (",perc*100,"%, n=",x[j,2],"/",nreads,")",sep=""),adj=0,cex=0.7)
					xlt=xlt-0.3
				}
			}
		}
	} else {
		text(5,xlt-0.4,"No DNA haplotype detected",adj=0,cex=0.7)
		xlt=xlt-0.3
	}
	xlt=xlt+0.15
	xl=min(xl,xlt)-1
	
	# writing rare variants on plot with gnomAD frequency
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
				text(0,xl,paste("Rare variant (AF<1e-5, VUS, hg19)*:",sep=""),adj=0,cex=1,col="orange")
				if(mut+WT>=10){	
					text(0,xl-0.3,paste(t1,", ",mut,"/",WT+mut," reads (",round(100*mut/(WT+mut),digits=1),"%) ","\n",t2,sep=""),adj=c(0,1),cex=0.8,col="orange")
				} else {
					text(0,xl-0.3,paste(t1,", ",mut,"/",WT+mut," reads (",round(100*mut/(WT+mut),digits=1),"%), warning: only ",WT+mut," reads","\n",t2,sep=""),adj=c(0,1),cex=0.8,col="orange")
				}
				xl=xl-1.4
			}
			if(pathoV==1){	
				
				if(NS[j,10]=="NM_020061.5:c.607T>C" | NS[j,10]=="NM_000513.2:c.607T>C"){
					Lmut=muts[which(muts[,1]==colnames(data3)[i] & muts[,2]=="C203R-LW"),3]
					Mmut=muts[which(muts[,1]==colnames(data3)[i] & muts[,2]=="C203R-MW"),3]
					Lwt=muts[which(muts[,1]==colnames(data3)[i] & muts[,2]=="C203-LW"),3]
					Mwt=muts[which(muts[,1]==colnames(data3)[i] & muts[,2]=="C203-MW"),3]
					text(0,xl,paste("Pathogenic variant (hg19)*:",sep=""),adj=0,cex=1,col="red")
					text(0,xl-0.3,paste(t1,", ",Lmut+Mmut,"/",Lmut+Mmut+Lwt+Mwt," reads (",round(100*(Lmut+Mmut)/(Lmut+Mmut+Lwt+Mwt),digits=1),"%) ","\n",t2,sep=""),adj=c(0,1),cex=0.8,col="red")
					text(0,xl-0.9,  paste(Lmut,"/",Lwt+Lmut," reads (",round(100*Lmut/(Lmut+Lwt+0.00001),digits=1),"%) on OPN1LW-exon4 and ",Mmut,"/",Mwt+Mmut," reads (",round(100*Mmut/(Mwt+Mmut+0.00001),digits=1),"%) on OPN1MW-exon4",sep="")  ,adj=c(0,1),cex=0.8,col="red")
					xl=xl-0.3
				} else {
					text(0,xl,paste("Pathogenic variant (hg19)*:",sep=""),adj=0,cex=1,col="red")
					text(0,xl-0.3,paste(t1,", ",mut,"/",WT+mut," reads (",round(100*mut/(WT+mut),digits=1),"%) ","\n",t2,sep=""),adj=c(0,1),cex=0.8,col="red")
				}
				xl=xl-1.4
			}
		}
	} else {
		text(0,xl,paste("No rare variants detected",sep=""),adj=0,cex=1,col=1)
		xl=xl-0.8
	}

	# writing deletions detected on the plot
	if(min(exons[,i])>0.1){
		text(0,xl,"No deletions detected",adj=0,col=1)
	} else {

		if(min(medi2)<10){
			text(3.2,xl,"(warning: some regions with very low median coverage)",adj=0,col="red",cex=0.7)
		}
		if(min(medi2)>=10 & min(medi2)<20){
			text(3.2,xl,"(warning: some regions with low median coverage)",adj=0,col="red",cex=0.7)
		}
		text(0,xl,"Deletion(s) detected:",adj=0,col="red",cex=1)

		xl=xl-0
		for (j in c(1:11)){
			if(exons[j,i]<0.1){
				text(0,xl-0.4,namesexons[j],adj=0,col="red",cex=0.8)
				xl=xl-0.4
			}
		}
	}

	# thr1=0.1 thr2=0.4 thr3=0.9

	# score1 = events above thr3 (1 each), can be spliceogenic haplo., patho. var., LCR del., full exon 3 or 6 del
	# score2 = events above thr2 (1 each), can be spliceogenic haplo. or patho. var.
	# score3 = events above thr1 (0.5 each), can be spliceogenic haplo. or patho. var.
	# score4 = rare VUS (score is percent of reads)
	# score5 dels L or M (1 per incomplete gene)
	# score6 Trp177Arg (if >90% of reads)
	# scoreA percentage of reads for spliceogenic haplotypes
	# scoreB percentage of reads for pathogenic variant
	
	# determination of phenotype based on events (score variables)
	pheno="Normal"
	if(score4>=thr1){pheno="Inconclusive"} # Rare VUS above 10%

	if(score2>=1 & score4>=thr1){pheno="Inconclusive (with extra VUS)"} # Event 40-90% with extra VUS above 10%
	if(score3>=1){pheno="Inconclusive"} # Event between 10 and 40%
	
	if(score3>=1 & score4>=thr1){pheno="Inconclusive (with extra VUS)"} # Event 10-40% with extra VUS above 10%

	if(score5>=1){pheno="Color vision deficiency suggested"} # Del L or M
	if(score5>=1 & score4>=thr1){pheno="Color vision deficiency suggested (with extra VUS)"}

	if(score2>=1 & scoreA>=thr2){pheno="Cone dysfunction disorder / BED suggested"} # Patho haplo 40-90%
	if(score2>=1 & scoreB>=thr2){pheno="Cone dysfunction disorder / BED suggested"} # Patho var 40-90%
	
	if(score5>=2){pheno="BCM suggested"} # del both L and M
	if(score2>=2){pheno="BCM suggested"} # Two events with 40-90% Susanne wrong!
	if(score1>=1){pheno="BCM suggested"} # One event >90%
	if(score1>=1 & score4>=thr1){pheno="BCM suggested (with extra VUS)"}
	if(score6>=1){pheno="Cone dystrophy suggested"} # Trp177Arg
	if(nLW==1 & nMW==1 & hybrid==0){pheno="BCM suggested"} # no complete gene

	# write phenotype on plot
	phenotype=rbind(phenotype,c(colnames(data3)[i],sex[i],score1+score2/2,pheno,round(LWs[i],digits=0),round(MWs[i],digits=0),CH,CH2+CH3))
	if(pheno=="Normal"){
		text(5,xl-1,paste("Inferred phenotype: ",pheno,sep=""),col=1,cex=1.1)
	}
	if(pheno=="Inconclusive (VUS)" | pheno=="Inconclusive"){
		text(5,xl-1,paste("Inferred phenotype: ",pheno,sep=""),col="orange",cex=1.1)
	}
	if(pheno!="Normal" & pheno!="Inconclusive (VUS)" & pheno!="Inconclusive"){
		text(5,xl-1,paste("Inferred phenotype: ",pheno,sep=""),col="red",cex=1.1)
	}

	## quality control
	# number of target having coverage between 0.1 and 0.7; it should not exist so indicate bad quality
	qual=length(which(ALL2[,i]>0.1 & ALL2[,i]<0.7))

	# deviation from copy numbers
	qual2=sum(abs(initial))
	# low quality warnings
	if(qual>2 & CH2==0){
		text(5,0,paste("Low quality, results not reliable!",sep=""),col="red",cex=1.2, xpd=NA)
		events=rbind(events,c(colnames(data3)[i],sex[i],"Low quality, results not reliable!","NA","NA","NA","NA","NA","NA"))
	}
	if(qual2>2.5 & (qual<=2 | CH2!=0) & sex[i]!="Female" & nLW>1 & nMW>1){
		text(5,0,paste("Copy-numbers not reliable",sep=""),col=1,cex=1.2, xpd=NA)
		events=rbind(events,c(colnames(data3)[i],sex[i],"Copy-numbers not fully reliable","NA","NA","NA","NA","NA","NA"))
	}
	if(sex[i]=="Female"){
		text(0,-0.55,paste("Warning: Hybrid gene prediction and copy-numbers in females is not fully reliable.",sep=""),adj=0,col=1,cex=0.7, xpd=NA)
	}
	if(sex[i]=="Undetermined"){
		text(0,-0.55,paste("Warning: Sex could not be predicted. Results can be unreliable.",sep=""),col="red",cex=0.7, xpd=NA,adj=0)
	}

	# second part of plotting with exon coverage
	beg=11.8
	m=max(2,ALL2*1,na.rm=T)
	m=max(2.5,ALL2[,i]*1.2,na.rm=T)
	m=max(m,5)
	par(mar=c(5, 4, 4, 2) + 0.1)
	plot(1:4,ylim=c(0,m),xlim=c(1,14.5), xaxt = "n",xlab="",ylab="Inferred copy-numbers",cex=0,main=paste("Copy-numbers of exons and deletions",sep=""),cex.main=1.2,cex.lab=1)
	for (j in 0:20){
		lines(c(-1,11.5),c(j,j),lty=2,col="darkgrey")
	}
	lines(c(-1,11.5),c(0,0))
	cols=c("#7c7c7d","#ff6f69","#ff6f69","#ff6f69","#ff6f69","#88d8b0","#88d8b0","#88d8b0","#88d8b0","grey","grey")
	for (j in 1:11){
		if(ALL2[j,i]<0.1){cols[j]="#8e07f5"}
	}
	for (j in 1:11){
		rect(j-0.3,0,j+0.3,ALL2[j,i],col=cols[j])
		if(cols[j]=="#8e07f5" & ALL2[j,i]<0.02){lines(c(j-0.3,j+0.3),c(ALL2[j,i],ALL2[j,i]),col="#8e07f5",lwd=6)}
	}
	n=dim(ALL2)[2]-1
	for (j in 1:11){
		selMale=which(sex=="Male" & colnames(ALL2)!=colnames(ALL2)[i])
		selFemale=which(sex=="Female" & colnames(ALL2)!=colnames(ALL2)[i])
		points(rep(j,length(selMale))+runif(length(selMale),-0.2,0.2),ALL2[j,selMale],pch=16,cex=0.4,col=1)
		points(rep(j,length(selFemale))+runif(length(selFemale),-0.2,0.2),ALL2[j,selFemale],pch=17,cex=0.4,col=4)
	}
	# adding information and legend on plot
	text(beg,m*0.98,"Inferred sex:",cex=1.2,adj=0)
	text(beg,m*0.93,sex[i],cex=1,adj=0)
	legend(beg,m*0.8,c("Deletion","Males","Females"),col=c("#8e07f5",1,4),pch=c(NA,16,17),lty=c(1,NA,NA),lwd=c(3,1,1),seg.len=1,pt.cex=c(1,1,1),bg="white",y.intersp=1.5,cex=1)
	text(11.8,m*0.05,"* LCR is usually not covered\n by WES and can have\n high variability.",cex=0.6,adj=0)
	# low quality regions
	low="No"
	par(xpd=NA)
	for (j in 1:11){
		if(medi2[j]<10){
			low="Yes"
			rect(j-0.3,m*(-0.075),j+0.3,m*(-0.17),border=NA,col="#ffcccc")
		}
		if(medi2[j]>=10 & medi2[j]<=20){
			low="Yes"
			rect(j-0.3,m*(-0.075),j+0.3,m*(-0.17),border=NA,col="#ffeead")
		}
	}
	axis(1, at=1:11, labels=rownames(ALL2),las=2,cex.axis=0.75)
	if(low=="Yes"){
		text(beg,m*(-0.07),"Exons with low coverage:",cex=0.7,adj=0)
		legend(beg,m*(-0.1),c("Very low cov. (<10)","Low cov. (<20)"),col=c("#ffcccc","#ffeead"),pch=c(15,15),seg.len=0.8,pt.cex=c(1,1),bg="white",y.intersp=1.5,cex=0.7)
	}
	lines(c(1.6,5.4),c(m*(-0.18),m*(-0.18)))
	lines(c(5.6,9.4),c(m*(-0.18),m*(-0.18)))
	lines(c(9.6,11.4),c(m*(-0.18),m*(-0.18)))
	text(3.5,m*(-0.21),"OPN1LW",cex=1,adj=0.5)
	text(7.5,m*(-0.21),"OPN1MW",cex=1,adj=0.5)
	text(10.5,m*(-0.21),"MW / LW",cex=1,adj=0.5)
	for (f in c(1.6,5.6,9.6,5.4,9.4,11.4)){lines(c(f,f),c(m*(-0.18),m*(-0.16)))}
	dev.off()

	# plot for probability of each copy number for LW and MW
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

# writing events and phenotypes to summary files
if(dim(events)[1]>0){
	for (i in 1:dim(events)[1]){
		events[i,3]=gsub("\n"," ",events[i,3])
	}
}
write.table(events,file=paste(dir,"/0.events.tsv",sep=""),quote=F,row.names=F,sep="\t")
write.table(phenotype[,c(1,2,4,5,6,7,8)],file=paste(dir,"/0.phenotypes.tsv",sep=""),quote=F,row.names=F,sep="\t")
