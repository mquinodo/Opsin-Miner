
#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

file1=args[1]
file2=args[2]
file3=args[3]

# read sequences
data=read.table(file=file1)

# convert into matrix
m=max(apply(data,1,nchar))
a=matrix(nrow=dim(data)[1],ncol=m)
for (i in 1:dim(data)[1]){
	temp=strsplit(as.character(data[i,1]),"")[[1]]
	a[i,1:length(temp)]=temp
}

ref=1:dim(a)[1]
for (i in 1:dim(a)[1]){
	sel=which(a[i,10:(dim(a)[2]-10)]=="G" & a[i,11:(dim(a)[2]-9)]=="C" & a[i,12:(dim(a)[2]-8)]=="T" & a[i,13:(dim(a)[2]-7)]=="G" & a[i,14:(dim(a)[2]-6)]=="C" & a[i,15:(dim(a)[2]-5)]=="A" & a[i,16:(dim(a)[2]-4)]=="T" & a[i,17:(dim(a)[2]-3)]=="C" & a[i,18:(dim(a)[2]-2)]=="A")+9
	if(length(sel)==1){
		ref[i]=sel
	} else {
		ref[i]=0
	}
}

# PROTEIN HAPLOTYPE
res=1:dim(a)[1]
for (i in 1:dim(a)[1]){
	if(ref[i]!=0 & ref[i]+9>0 & ref[i]+30<=length(a[i,]) & is.na(a[i,ref[i]+26])==F & dim(a)[2]>=ref[i]+26){
		# 230
		if(a[i,ref[i]+9]=="T"){
			res[i]="I"
		} else if(a[i,ref[i]+9]=="C"){
			res[i]="T"
		} else {
			res[i]="X"
		}
		# 233
		if(a[i,ref[i]+17]=="G" & a[i,ref[i]+18]=="C" & a[i,ref[i]+19]=="T"){
			res[i]=paste(res[i],"A",sep="")
		} else if(a[i,ref[i]+17]=="A" & a[i,ref[i]+18]=="G" & a[i,ref[i]+19]=="C"){
			res[i]=paste(res[i],"S",sep="")
		} else {
			res[i]=paste(res[i],"X",sep="")
		}
		# 236
		if(a[i,ref[i]+26]=="A"){
			res[i]=paste(res[i],"M",sep="")
		} else if(a[i,ref[i]+26]=="G"){
			res[i]=paste(res[i],"V",sep="")
		} else {
			res[i]=paste(res[i],"X",sep="")
		}
	} else {
		res[i]="NA"
	}
}
haplotype=res[which(res!="NA")]

out=unique(haplotype)
out=cbind(out,out,out,out)
if(dim(out)[1]>0){for (i in 1:dim(out)[1]){
	out[i,2]=length(which(haplotype==out[i,1]))
	out[i,3]=length(haplotype==out[i,1])
	out[i,4]=round(as.numeric(out[i,2])/as.numeric(out[i,3]),3)
}}
colnames(out)=c("Haplotype","Reads","Total","Ratio")

write.table(out,file=file2,quote=F,row.names=F,sep="\t")
