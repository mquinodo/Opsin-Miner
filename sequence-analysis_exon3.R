
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
	sel=which(a[i,10:(dim(a)[2]-10)]=="A" & a[i,11:(dim(a)[2]-9)]=="G" & a[i,12:(dim(a)[2]-8)]=="A" & a[i,13:(dim(a)[2]-7)]=="T" & a[i,14:(dim(a)[2]-6)]=="T" & a[i,15:(dim(a)[2]-5)]=="T" & a[i,16:(dim(a)[2]-4)]=="G" & a[i,17:(dim(a)[2]-3)]=="A" & a[i,18:(dim(a)[2]-2)]=="T")+9
	if(length(sel)==1){
		ref[i]=sel
	} else {
		ref[i]=0
	}
}

# PROTEIN HAPLOTYPE
res=1:dim(a)[1]
for (i in 1:dim(a)[1]){
	if(ref[i]!=0 & ref[i]-30>0 & ref[i]+51<=length(a[i,])){
		if(is.na(a[i,ref[i]+51])==F){
			# 153
			if(a[i,ref[i]-30]=="A"){
				res[i]="M"
			}
			if(a[i,ref[i]-30]=="C"){
				res[i]="L"
			}
			if(res[i]==i){
				res[i]="X"
			}
			# 171
			if(a[i,ref[i]+24]=="G"){
				res[i]=paste(res[i],"V",sep="")
			}
			if(a[i,ref[i]+24]=="A"){
				res[i]=paste(res[i],"I",sep="")
			}
			if(a[i,ref[i]+24]!="G" & a[i,ref[i]+24]!="A"){
				res[i]=paste(res[i],"X",sep="")
			}
			# 174
			if(a[i,ref[i]+34]=="T"){
				res[i]=paste(res[i],"V",sep="")
			}
			if(a[i,ref[i]+34]=="C"){
				res[i]=paste(res[i],"A",sep="")
			}
			if(a[i,ref[i]+34]!="T" & a[i,ref[i]+34]!="C"){
				res[i]=paste(res[i],"X",sep="")
			}
			# 178
			if(a[i,ref[i]+45]=="G"){
				res[i]=paste(res[i],"V",sep="")
			}
			if(a[i,ref[i]+45]=="A"){
				res[i]=paste(res[i],"I",sep="")
			}
			if(a[i,ref[i]+45]!="G" & a[i,ref[i]+45]!="A"){
				res[i]=paste(res[i],"X",sep="")
			}
			# 180
			if(a[i,ref[i]+51]=="G"){
				res[i]=paste(res[i],"A",sep="")
			}
			if(a[i,ref[i]+51]=="T"){
				res[i]=paste(res[i],"S",sep="")
			}
			if(a[i,ref[i]+51]!="G" & a[i,ref[i]+51]!="T"){
				res[i]=paste(res[i],"X",sep="")
			}
		} else {
		res[i]="NA"
		}
	} else {
		res[i]="NA"
	}
}
haplotype=res[which(res!="NA")]

# DNA HAPLOTYPE
res=1:dim(a)[1]
for (i in 1:dim(a)[1]){
	if(ref[i]!=0 & ref[i]-34>0 & ref[i]+51<=length(a[i,])){
		if(is.na(a[i,ref[i]+51])==F){
			# 151
			res[i]=a[i,ref[i]-34]
			# 153
			res[i]=paste(res[i],a[i,ref[i]-30],sep="")
			# 155
			res[i]=paste(res[i],a[i,ref[i]-22],sep="")
			# 171
			res[i]=paste(res[i],a[i,ref[i]+24],sep="")
			# 171-2
			res[i]=paste(res[i],a[i,ref[i]+26],sep="")
			# 174
			res[i]=paste(res[i],a[i,ref[i]+34],sep="")
			# 178
			res[i]=paste(res[i],a[i,ref[i]+45],sep="")
			# 180
			res[i]=paste(res[i],a[i,ref[i]+51],sep="")
		} else {
			res[i]="NA"
		}
	} else {
		res[i]="NA"
	}
}
haplotypeDNA=res[which(res!="NA")]

# rs949430 rs713 rs731614 rs5986963 rs5986964 rs149897670 rs145009674 rs949431
# R151R L153M V155V V171M V171V A174V I178V S180A

out=unique(haplotype)
out=cbind(out,out,out,out)
for (i in 1:dim(out)[1]){
	out[i,2]=length(which(haplotype==out[i,1]))
	out[i,3]=length(haplotype==out[i,1])
	out[i,4]=round(as.numeric(out[i,2])/as.numeric(out[i,3]),3)
}
colnames(out)=c("Haplotype","Reads","Total","Ratio")

write.table(out,file=file2,quote=F,row.names=F,sep="\t")

out=unique(haplotypeDNA)
out=cbind(out,out,out,out)
if(dim(out)[1]>0){for (i in 1:dim(out)[1]){
	out[i,2]=length(which(haplotypeDNA==out[i,1]))
	out[i,3]=length(haplotypeDNA==out[i,1])
	out[i,4]=round(as.numeric(out[i,2])/as.numeric(out[i,3]),3)
}}
colnames(out)=c("Haplotype","Reads","Total","Ratio")

write.table(out,file=file3,quote=F,row.names=F,sep="\t")
