#writed by tongxinkai

## args[1] pwd
## args[2] gene id

rm(list=ls());options(stringsAsFactors=F)
args <- commandArgs(trailingOnly = TRUE)
setwd(paste(args[1],args[2],sep="/"))
#a=read.table('chr14-27853347.ld',header=F)
#colnames(a)=c('ps','r2')
b=read.table('ttt_7',header=F)
h=read.table('cs7.bed',header=F)
colnames(h)=c('ps','sc')
h$ps=h$ps/1e6

colnames(b)=c('gene','chr','ps','p')
#dat=merge(a,b,by='ps')
dat=b
dat$ps=dat$ps/1e6;dat$logp=-log(dat$p,10)
dat1=subset(dat,logp<=100)
dat2=subset(dat,logp>100)
dat2=dat1[1,]
#dat3=subset(dat,r2>0.4 & r2<=0.6)
#dat4=subset(dat,r2>0.6 & r2<=0.8)
#dat5=subset(dat,r2>0.8 & r2<1)
#dat6=subset(dat,r2==1)
ymin=min(dat$logp);ymax=max(dat$logp)
gene=read.csv("g.csv",header=T)
gene[,c(2,3)]=gene[,c(2,3)]/1e6
for(i in 1:nrow(gene)){
  if(gene[i,4]==''){
    n=nchar(gene[i,1])
    gene[i,4]=substr(gene[i,1],n-5,n)
  }
}
gene=gene[,-1]
colnames(gene)=c('st','end','name');gene=gene[order(gene$st),]
source("../gene.r")
l=delovlp(gene);n=length(l)
m1=ifelse(min(dat$ps)>min(gene$st),min(gene$st),min(dat$ps))
m2=ifelse(max(dat$ps)>max(gene$end),max(dat$ps),max(gene$end))
png(paste(args[1],"/",args[2],"/",args[2],"_tmp77.png",sep=""),height=6000,width=6000,res=700)
#plot(dat1$logp~dat1$ps,xlab='position',ylab='-log(p)',cex.lab=1.3,xlim=range(dat$ps),
     #ylim=c(0,ymax*1.2),xaxt='n',col='navy',cex=0.8)

#par(fig=c(0,1,0,0.7))
layout(matrix(c(1,2),2,1,byrow=TRUE),widths=c(3,3),heights=c(1,2))	 	 

plot(h$sc~h$ps,xlab="",ylab='CAVIAR score',main=paste("STR ",dat1[1,2],"_",dat1[1,3]*1e6," (",dat1[1,4],")",sep=""),cex.lab=1.5,cex.main=1.5,col='navy',cex=0.7,ylim=c(0,1),xlim=c(m1-0.01,m2+0.01))
points(h[1,1],h[1,2],col='red',cex=1.25,pch=15)


plot(dat1$logp~dat1$ps,xlab='Position (Mb)',ylab='-log (pvalue)',main="",cex.lab=1.5,
     ylim=c(-5-5*n,ymax*1.2),
     xlim=c(m1-0.01,m2+0.01),col='navy',cex=0.7)
#axis(side=2,at=c(-4,-2,0,2,4,6,8),labels=c(-4,-2,0,2,4,6,8), cex.lab=1.3)
abline(h=0)
#points(dat2$logp~dat2$ps,col='red',cex=0.75,pch=15)
points(dat1[1,3],dat1[1,5],col='red',cex=1.25,pch=15)
legend('bottomleft',col=c("red","blue"),pch=c(0,1),legend=c("STR","SNP"),cex=1.15)
#points(dat3$logp~dat3$ps,col='olivedrab2',cex=0.5)
#points(dat4$logp~dat4$ps,col='orange',cex=0.5)
#points(dat5$logp~dat5$ps,col='red3',cex=0.5)
#points(dat6$logp~dat6$ps,col='red3',cex=1)
#text(dat6$ps+0.04,dat6$logp+0.6,labels = 'rs341574911',cex=0.8)
#text(dat6$ps+0.14,dat6$logp+0.3,labels = 'rs330520570',cex=0.6)
#text(dat6$ps+0.14,dat6$logp,labels = 'rs324538001',cex=0.6)
#abline(h=0)
#axis(1,at=seq(24.8,25.9,0.2),labels=seq(24.8,25.9,0.2),3)
#legend('topleft',pch=1,col=c('red3','orange','olivedrab2','steelblue1','navy'),
#       legend=c('r2>0.8','0.6<r2<0.8','0.4<r2<0.6','0.2<r2<0.4','0.2<r2'),cex=0.65)
for(i in 1:n){
  rect(xleft = l[[i]]$st,
       ybottom = rep(0-i*5,length(l[[i]]$st)),
       xright =l[[i]]$end,
       ytop=rep(1-i*5,length(l[[i]]$st)),
       col='green',border = 'green')
  text((l[[i]]$st+l[[i]]$end)/2,-2-i*5,labels = l[[i]]$name,cex=0.8)
}

#par(fig=c(0,1,0.7,1))

#plot(dat1$sc~dat1$ps,xlab="",ylab='CAVIAR score',main=paste("STR ",dat1[1,2],"_",dat1[1,3]*1e6," (",dat1[1,4],")",sep=""),cex.lab=1.5,cex.main=1.5,col='navy',cex=0.7,ylim=c(0,1),xlim=c(m1-0.01,m2+0.01))
#points(dat1[1,3],dat1[1,5],col='red',cex=1.25,pch=15)

dev.off()

pdf(paste(args[1],"/",args[2],"/",args[2],"_tmp77.pdf",sep=""),height=9,width=9)
#plot(dat1$logp~dat1$ps,xlab='position',ylab='-log(p)',cex.lab=1.3,xlim=range(dat$ps),
     #ylim=c(0,ymax*1.2),xaxt='n',col='navy',cex=0.8)

#par(fig=c(0,1,0,0.7))
layout(matrix(c(1,2),2,1,byrow=TRUE),widths=c(3,3),heights=c(1,2))

plot(h$sc~h$ps,xlab="",ylab='CAVIAR score',main=paste("STR ",dat1[1,2],"_",dat1[1,3]*1e6," (",dat1[1,4],")",sep=""),cex.lab=1.5,cex.main=1.5,col='navy',cex=0.7,ylim=c(0,1),xlim=c(m1-0.01,m2+0.01))
points(h[1,1],h[1,2],col='red',cex=1.25,pch=15)

	 	 
plot(dat1$logp~dat1$ps,xlab='Position (Mb)',ylab='-log (pvalue)',main="",cex.lab=1.5,
     ylim=c(-5-5*n,ymax*1.2),
     xlim=c(m1-0.01,m2+0.01),col='navy',cex=0.7)
#axis(side=2,at=c(-4,-2,0,2,4,6,8),labels=c(-4,-2,0,2,4,6,8), cex.lab=1.3)
abline(h=0)
#points(dat2$logp~dat2$ps,col='red',cex=0.75,pch=15)
points(dat1[1,3],dat1[1,5],col='red',cex=1.25,pch=15)
legend('bottomleft',col=c("red","blue"),pch=c(0,1),legend=c("STR","SNP"),cex=1.15)
#points(dat3$logp~dat3$ps,col='olivedrab2',cex=0.5)
#points(dat4$logp~dat4$ps,col='orange',cex=0.5)
#points(dat5$logp~dat5$ps,col='red3',cex=0.5)
#points(dat6$logp~dat6$ps,col='red3',cex=1)
#text(dat6$ps+0.04,dat6$logp+0.6,labels = 'rs341574911',cex=0.8)
#text(dat6$ps+0.14,dat6$logp+0.3,labels = 'rs330520570',cex=0.6)
#text(dat6$ps+0.14,dat6$logp,labels = 'rs324538001',cex=0.6)
#abline(h=0)
#axis(1,at=seq(24.8,25.9,0.2),labels=seq(24.8,25.9,0.2),3)
#legend('topleft',pch=1,col=c('red3','orange','olivedrab2','steelblue1','navy'),
#       legend=c('r2>0.8','0.6<r2<0.8','0.4<r2<0.6','0.2<r2<0.4','0.2<r2'),cex=0.65)
for(i in 1:n){
  rect(xleft = l[[i]]$st,
       ybottom = rep(0-i*5,length(l[[i]]$st)),
       xright =l[[i]]$end,
       ytop=rep(1-i*5,length(l[[i]]$st)),
       col='green',border = 'green')
  text((l[[i]]$st+l[[i]]$end)/2,-2-i*5,labels = l[[i]]$name,cex=0.8)
}

#par(fig=c(0,1,0.7,1))

#points(dat1$sc~dat1$ps,xlab='',ylab='CAVIAR score',main=paste("STR ",dat1[1,2],"_",dat1[1,3]*1e6," (",dat1[1,4],")",sep=""),cex.lab=1.5,cex.main=1.5,col='navy',cex=0.7,ylim=c(0,1),xlim=c(m1-0.01,m2+0.01))
#points(dat1[1,3],dat1[1,5],col='red',cex=1.25,pch=15)

dev.off()

