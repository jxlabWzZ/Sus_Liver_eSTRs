 ### shell script to run peer
 
 Rscript=/home/wuzhongzi/miniconda3/envs/py27/bin/Rscript
 for j in {5,10,15,20,25,30,35,40,45,50,55};do
	 sed "s/PEER_setNk(model,15)/PEER_setNk(model,$j)/g" tmp.r > tmp$j.r
	 sed "s/i=1/i=$j/g" tmp$j.r > round$j.r
	 nohup ${Rscript} round$j.r &
	 done

	 
	 
###
### tmp.r Rscript

##/home/wuzhongzi/miniconda3/envs/py27/bin/R
rm(list=ls());options(stringsAsFactors=FALSE)
setwd("/home/wuzhongzi/hip_str/RNAdata/14.rpkm")
load("RNASeq_TPM.RData")
data=RNASeq_TPM
library(peer)
i=1
model = PEER()
PEER_setPhenoMean(model,as.matrix(data))
PEER_setNk(model,15)#infer K=15 hidden confounders,
PEER_update(model)
factors = PEER_getX(model)
save(factors,file=paste("./","RNASeq_factor_",i,".Rdata",sep=""))
write.table(factors, file=paste("./","RNASeq_factors_",i,".txt",sep=""), sep="\t")



a<-read.table("../../F6F7_556_Liver_FPKM.bed",header=T)
b<-a[,-1]
RNASeq_TPM<-t(b)
save.image("RNASeq_F67_TPM.Rdata")



load("RNASeq_F6_TPM.RData")
d<-rownames(RNASeq_TPM)

load("RNASeq_factor_45.Rdata")
rownames(factors)<-d
write.table(factors,file="f6_p45.txt",sep="\t",row.names=T,col.names=F,quote=FALSE)


load("RNASeq_F7_TPM.RData")
d<-rownames(RNASeq_TPM)

load("RNASeq_factor_45.Rdata")
rownames(factors)<-d
write.table(factors,file="f7_p45.txt",sep="\t",row.names=T,col.names=F,quote=FALSE)


for i in 1 5 10 15 20 25 30 35 40 45
do
less -S f6_p${i}.txt |sed -n '1,1p' | awk 'BEGIN{OFS="\t"}{$1="ID"; for(i=2;i<=NF;i++) $i="PE"(i-1); print $0}' | cat - f6_p${i}.txt | paste f6_260_cov0_01.txt - | cut -f 1-18,20- > cov_f${i}.txt
done



for i in 1 5 10 15 20 25 30 35 40 45
do
less -S f7_p${i}.txt |sed -n '1,1p' | awk 'BEGIN{OFS="\t"}{$1="ID"; for(i=2;i<=NF;i++) $i="PE"(i-1); print $0}' | cat - f7_p${i}.txt | paste f7_296_cov0_01.txt - | cut -f 1-19,21- > cov_f${i}.txt
done



