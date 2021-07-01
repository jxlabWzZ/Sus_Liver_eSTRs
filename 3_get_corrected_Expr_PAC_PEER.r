## /home/wuzhongzi/miniconda3/envs/py27/bin/R
## conda activate py27 (r-3.4.1)
## conda install -c bioconda bioconductor-preprocesscore

rm(list=ls());options(stringsAsFactors=F)
library(preprocessCore)

Normalize <- function(m){ apply(m,2,function(x){(x - mean(x))/sd(x)} ) } #Using normalize.quantiles.use.target instead
##	Function that generate the Hat Matrix
HatMatrix <- function(c){ return (c %*% solve(t(c) %*% c, tol=1e-23) %*% t(c))}
##	Function that generates the residuals matrix
GenerateResidual <- function(h,y){ return((diag(nrow(h)) - h) %*% y)}


# Building covariate matrix

##	 Open covariate files & set up variables
# pcas <- read.table("/storage/szfeupe/Runs/GTEx_estr/gtex.pca", sep=" ", header=FALSE)**
# pcas include batch, age, weight and sex et. al.
# peer factors 0,1,3,5,7,10,15,20,25,30,35,40,45,50,55 factors
#

Covar <- read.table("/home/wuzhongzi/hip_str/RNAdata/15.eSTR/cov_f30.txt", sep="\t", header=T)

head(Covar[,c(1:6)])

# id   s4357   s4364   s4365   s4372   s4375
# 1  Gen   0.000   0.000   0.000   0.000   0.000
# 2  Sex   0.000   1.000   0.000   1.000   0.000
# 3 Data 255.000 255.000 255.000 255.000 255.000
# 4   Wt  68.638  76.057  54.144  81.233  78.786

## 	Transform covariate to matrix and call it C, (SxC) 
C<-data.matrix(t(Covar[,-1]))

#	Load Expression data
#	Check header and rownames 
#	Basically you want rownames(Expr) to be the ENSGxxxx ids and colnames(Expr) to be the Ids of samples GTEXyyyyy
#	if not change it like so [eg. rownames(Expr) <- Expr$GeneIDs and then Expr$GeneIDs = NULL ]

b<- read.table("/home/wuzhongzi/hip_str/RNAdata/15.eSTR/rpkm548.txt", sep="\t", header=T)
a<-read.table("/home/wuzhongzi/hip_str/RNAdata/15.eSTR/rpkm548.gene3", sep="\t", header=F)
Expr<-b[,-1]
rownames(Expr)<-t(a)
dim(Expr)
Yo=Expr[,rownames(C)]

##	Transpose expression matrix to get samples in rows (SxC)
Y <- data.matrix(t(Yo))
dim(Y)

##	CORRECTION
#	Let's get the Hat matrix
H <- HatMatrix(C)
dim(H)

#Fitting expression into normal dist mean=0, sd=1
Y1=normalize.quantiles.use.target(as.matrix(Y), rnorm(nrow(Y)))
dimnames(Y1) = dimnames(Y)
head(Y1[,1])

## 	Let's get the residual expression matrix Y.prime 
Y.prime <- GenerateResidual(H,Y1)
dim(Y.prime)
#head(rownames(Y.prime))
colnames(Y.prime)<-rownames(Expr)
##	Record it
write.table(Y.prime, file="Corr_Expr_f30.csv", sep=",")

## write.table(Y.prime, file="Corr_Expr.txt",row.names=T,col.names=T,quote = F, sep="\t")

# We can ADD population covariates, here: age and sex but there are many missing data
#
#
q()

