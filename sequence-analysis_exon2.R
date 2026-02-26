
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
	sel=which(a[i,10:(dim(a)[2]-10)]=="T" & a[i,11:(dim(a)[2]-9)]=="G" & a[i,12:(dim(a)[2]-8)]=="C" & a[i,13:(dim(a)[2]-7)]=="A" & a[i,14:(dim(a)[2]-6)]=="T" & a[i,15:(dim(a)[2]-5)]=="C" & a[i,16:(dim(a)[2]-4)]=="C" & a[i,17:(dim(a)[2]-3)]=="G" & a[i,18:(dim(a)[2]-2)]=="T")+9
	if(length(sel)==1){
		ref[i]=sel
	} else {
		ref[i]=0
	}
}

# PROTEIN HAPLOTYPE
res=1:dim(a)[1]
for (i in 1:dim(a)[1]){
	if(dim(a)[2]>=ref[i]+152){
		if(ref[i]!=0 & ref[i]+2>0 & ref[i]+30<=length(a[i,]) & is.na(a[i,ref[i]+152])==F){
			# 65
			if(a[i,ref[i]-1]=="C"){
				res[i]="T"
			} else if(a[i,ref[i]-1]=="T"){
				res[i]="I"
			} else {
				res[i]="X"
			}
			# 110
			if(a[i,ref[i]+136]=="A"){
				res[i]=paste(res[i],"I",sep="")
			} else if(a[i,ref[i]+136]=="G"){
				res[i]=paste(res[i],"V",sep="")
			} else {
				res[i]=paste(res[i],"X",sep="")
			}
			# 116
			if(a[i,ref[i]+152]=="C"){
				res[i]=paste(res[i],"S",sep="")
			} else if(a[i,ref[i]+152]=="A"){
				res[i]=paste(res[i],"Y",sep="")
			} else {
				res[i]=paste(res[i],"X",sep="")
			}
		} else {
			res[i]="NA"
		}
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
