library(MASS)
library(TExPosition)
library(ExPosition)
library(ggplot2)                #change zero position in data simulation to reflect orthogonality between components. 
library(reshape2)
source("SimulationGen/catConvert.R")
source("PLSR/GPLSos.R")
source("VariableSel_PLSCA/InternalFunc/plscaGPLS.R")
source("VariableSel_PLSCA/InternalFunc/plscaGPLSStep1.R")
source("VariableSel_PLSCA/InternalFunc/plscaGPLSStep2.R")
source("VariableSel_PLSCA/InternalFunc/rowNormsvs.R")
source("VariableSel_PLSCA/InternalFunc/caSVD.R")
source("VariableSel_PLSCA/InternalFunc/makeRowProfilesvs.R")
source("VariableSel_PLSCA/InternalFunc/caNormvs.R")
source("CatPLS/createIndMatrix.R")
library(sgPLS)
library(mvtnorm)
library(Matrix)
###########Define Matrix Size############plevel
n <- 100                        #Define size of observation
p <- 20                        #No. of X variables
q <- 3                         #No. of Z variables
off.var <- 5                   #No. of variables not relevant
ncomp <- 2
set.seed(36)                   #determine the levels in each cat variable randomly
plevel <- rep(3,p)
qlevel <- rep(3,q)
#plevel <- sample(seq(3,4,by=1),p,replace=TRUE)  #Generate no. of levels for X variable
#qlevel <- sample(seq(3,4,by=1),q,replace=TRUE)  #Generate no. of levels for X variable
###########End Size Defi############
###########Built Data Matrix##############
data.raw.com <- matrix(rnorm(n*ncomp,0,1),n,ncomp)     #generate the matrix with random variables from normdis with size n-off.var by n
#data.raw.c <- matrix(rnorm(ncomp*(p-off.var),0,1),ncomp,(p-off.var))                         #for generate x data
#data.raw.c <- cbind(data.raw.c,c(rep(0,ncomp)),c(rep(0,ncomp)))      #Add two 0 rows in the matrix to ensure the correlation between X and Z variables does not relate to the last to variables 
#data.raw.d <- matrix(rnorm(ncomp*q,0,1),ncomp,q)                         #for generate z data
data.raw.c <- matrix(c(sample(c(-1,1,1.5),(p-off.var),replace=TRUE),rep(0,off.var)),ncomp,p,byrow=TRUE)
data.raw.d <- matrix(runif(ncomp*q,0.5,2),ncomp,q)

rho=0.5
sigmax=matrix(0,p-off.var,p-off.var)
sigmax.noise=matrix(0,off.var,off.var)
diag(sigmax.noise)<-1
for(i in 1:(p-off.var)){
  for(j in 1:(p-off.var)){
    sigmax[i,j]<-rho^(abs(i-j))
  }
}
sigma.x=bdiag(sigmax,sigmax.noise)
sigma.x=sigma.x*(rho**2)
sigma.x=as.matrix(sigma.x)

sigmay=matrix(0,q,q)
for(i in 1:q){
  for(j in 1:q){
    sigmay[i,j]<-rho^(abs(i-j))
  }
}
sigma.y=sigmay*(rho**2)
sigma.y=as.matrix(sigma.y)

data.raw.errx <- rmvnorm(n,mean=rep(0,p),sigma=sigma.x,method="svd")
data.raw.erry <- rmvnorm(n,mean=rep(0,q),sigma=sigma.y,method="svd")



#data.raw.erc <- matrix(sample(c(-0.1,0.1),n*off.var,replace=TRUE),n,off.var)                       #error for x data
#data.raw.erd <- matrix(sample(c(-0.15,0.13),n*q,replace=TRUE),n,q)                       #error for z data
data.raw.x <- data.raw.com %*% data.raw.c + (data.raw.errx)
data.raw.z <- data.raw.com %*% data.raw.d + (data.raw.erry)
#data.raw.x[,(p-off.var+1):p] <- data.raw.x[,(p-off.var+1):p]+ (data.raw.erc*0.01)

x.cat <- catConvert(data.raw.x,plevel)
z.cat <- catConvert(data.raw.z,qlevel)
###########End Build######################
###########Define Result Structure########
#res.os <- osPlsr(x.cat,z.cat,2,1/sqrt(n))

data1=createIndMatrix(x.cat,off.var,plevel)
data2=createIndMatrix(z.cat,0,qlevel)


#create ind.block.x and ind.block.z
ind.block.x<-ind.block.z<-NULL
ind.block.x<-plevel[1]
ind.block.z<-qlevel[1]
for(i in 1:(p-off.var)){
  ind.block.x=cbind(ind.block.x,sum(plevel[1:(i+1)]))
}
for(i in 1:(q-2)){
  ind.block.z=cbind(ind.block.z,sum(qlevel[1:(i+1)]))
}
data1[,(ind.block.x[p-2]+1):dim(data1)[2]]=data1[,(ind.block.x[p-2]+1):dim(data1)[2]]
#end create ind.block.x and z
#Mass and weight
R.ini <- t(data1)%*%data2
colTotal <- colSums(R.ini)
rowTotal <- rowSums(R.ini)
rowTotal[which(rowTotal==0)]<-1
grandTotal <- sum(R.ini)

masses = rowTotal/grandTotal
#masses<-ifelse(masses==0,0.001,masses)
weights = (colTotal/grandTotal)^(-1)


#res.ca <- plscaConv(DATA1=data1,DATA2=data2,k=3)
res.ca.vs <- plscaGPLS(data1, data2, ncomp = 2, mode = "regression", keepX = (rep(p-off.var,p)), 
                   keepZ = (rep(q,q)), ind.block.x = ind.block.x
                   ,ind.block.z = ind.block.z, masses = masses, weights=weights)

###########END Result Structure##########
###########Plot Res Comparison###########
# plotdata<-data.frame(
#   Xi = res.os$Xi[,1],
#   lx = res.ca$TExPosition.Data$lx[,1],
#   xvalue = seq(from=1,to=nrow(data.raw.x),by=1))
# plotdatalong <- melt(plotdata,id="xvalue")
# lvxresults <- ggplot(plotdatalong,
#                   aes(x=xvalue,y=value,colour=variable,shape=variable))+geom_point()+geom_line()+scale_shape_manual(values=c(1,2))+labs(title="Scale changed",x="",y="") 
# print(lvxresults)
# plotdata<-data.frame(
#   Omega = res.os$Omega[,1],
#   ly = res.ca$TExPosition.Data$ly[,1],
#   xvalue = seq(from=1,to=nrow(data.raw.x),by=1))
# plotdatalong <- melt(plotdata,id="xvalue")
# lvyresult <- ggplot(plotdatalong,
#                   aes(x=xvalue,y=value,colour=variable,shape=variable))+geom_point()+geom_line()+scale_shape_manual(values=c(1,2))+labs(title="Scale changed",x="",y="") 
# print(lvyresult)
###########Plot Res Comparison###########
