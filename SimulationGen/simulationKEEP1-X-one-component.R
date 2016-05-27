###################################################
### Simu 1: as Lin (2013) paper but 1 component + sequence noise
###################################################
#source("Dropbox/PROJECT-A-STAR-PLS/sPLS.BP.R")
#source("Dropbox/PROJECT-A-STAR-PLS/perf-BP.R")
require(Matrix)
require(mvtnorm)
require(xtable)
folds=5
n <- 100
sigma.gamma <- 1
p <- 400
q <- 10
ind.block.x <- c(0,seq(20,380,20))                # each block has 20 variables
set.seed(5)
j <- sample(ind.block.x,4)                        # select the start indices for non-zero blocks
position.block.x <-j/20+1                         # as above
theta.xx <- rep(0,400)                      
theta.x <- c(rep(1,15),rep(-1,30),rep(1.5,15))    # only 15 non-zero entries for SGPLS application
k <- 1
for (ii in 1:4){
  theta.xx[(j[ii]+1):(j[ii]+15)] <- theta.x[(k):(k+14)] #Insert non-zero theta.x into theta.xx at predefined indices
  k <- k+15
}
theta.xx1 <- theta.xx

folds = 5
set.seed(10)
theta.yy1 <- runif(q,0.5,2)
#theta.yy2 <- runif(q,0.5,2)
cov.irrelevant <- 0#runif(1,0,0.2)
rho <- 0.5
ng.x <- 20                            #no.of group of x
nx <- 20                              #no.of x in each group
set.seed(29)                          
rho1 <- runif(1,0.2,0.4)              #
rho1 <- 0.5
set.sigma.e <- c(0.5,2.5) #seq(1.5,2.5,length=5)
res.qua.group <- res.qua.spls <- res.qua.sparse <- NULL
for (kk  in 1:2){
  sigma.e <- set.sigma.e[kk]
#  sigma.e <- 0.5
#  print(c(sigma.e,kk))
  
  Sigma.noise <- matrix(cov.irrelevant/(sigma.e**2),nrow=5,ncol=5)
  diag(Sigma.noise) <- 1
  
  rho <- 0.5
  ng.x <- 20
  nx <- 20
  Sigma <- matrix(cov.irrelevant/sigma.e,nrow=15,ncol=15)
  for (i in 1:(ng.x-5)) {
    for(j in 1:(ng.x-5)){
      Sigma[i,j] <- rho^(abs(i-j))
    }
  }
  diag(Sigma) <- 1
  
  Sigma <- bdiag(Sigma,Sigma.noise)          #Combine the diagnal of both matrix to form a new matrix
  
  
  L <- list()
  for (i in position.block.x){              #using sigma to form a list, current error is generated to use for non-zero x groups
    L[[i]] <- Sigma 
  }
  
  Sigma.noise <- matrix(cov.irrelevant/(sigma.e**2),nrow=20,ncol=20)
  diag(Sigma.noise) <- 1
  
  for (i in setdiff(1:20,position.block.x)){     #set error to be an identity matrix for other block
    L[[i]] <- Sigma.noise 
  }
  
  
  Sigma.x <- do.call(bdiag,L)                    #using each list within L diagnoally to form a matrix
  Sigma.x[which(Sigma.x==0,arr.ind=TRUE)] <- 0
  Sigma.x <- Sigma.x*(sigma.e**2)
  Sigmax <- as.matrix(Sigma.x)

rho <- 0.5
Sigma <- matrix(nrow=q,ncol=q)
for (i in 1:q) {
  for(j in 1:q){
    Sigma[i,j] <- rho^(abs(i-j))
  }
}


Sigma.y <- Sigma*(sigma.e**2)
Sigmay <- as.matrix(Sigma.y)

ind.block.x <- seq(20,380,20)

trueX <- which(theta.xx1!=0)

res.groupX1 <- res.sparseX1 <- res.splsX1 <- NULL
#res.groupX2 <- res.sparseX2 <- res.splsX2 <- NULL

grid.X <- c(40,50,60,70,80,100,140,180,250,300)
#grid.X <- c(40,50)
grid.gX <- c(2,3,4,5,6,8,10,15,20)
grid.gX <- c(4)
grid.alpha.X <- c(0.05,0.25,0.5,0.75,0.95)
grid.alpha.X <- seq(0.05,0.95,length=20)
#grid.alpha.X <- c(0.05,0.25)
choicesetseed <- 15
Nsimu <- 1
tab.group <- tab.sparse <- tab.spls <- matrix(ncol=3,nrow=Nsimu)
for (isimu in 1:Nsimu) {
  gam1 <- rnorm(n,0,1)
#  gam2 <- rnorm(n,0,1)
  print(isimu)
  X <- matrix(c(gam1),ncol=1,byrow=FALSE)%*%matrix(c(theta.xx1),nrow=1,byrow=TRUE)+rmvnorm(n,mean=rep(0,p),sigma=Sigmax,method="svd")
  Y <- matrix(c(gam1),ncol=1,byrow=FALSE)%*%matrix(c(theta.yy1),nrow=1,byrow=TRUE)+rmvnorm(n,mean=rep(0,q),sigma=Sigmay,method="svd")
  
####Group PLS
#res.tun.gpls <- choice.tuning.gpls.X.comp(X,Y,folds = folds,ncomp=1,grid.gX,setseed=choicesetseed,ind.block.x=ind.block.x)
model.group <- Group.spls.BP(X,Y,ncomp=1,mode="regression",max.iter=500,tol=1e-06,keepX=rep(4,1),ind.block.x=ind.block.x)

if(norm(model.group$loadings$X[,1]-theta.xx1) < norm(-model.group$loadings$X[,1]-theta.xx1)) temp <- model.group$loadings$X[,1] else  temp <- -model.group$loadings$X[,1]
res.groupX1 <- cbind(res.groupX1,temp)


  
quality <- quality.simu.X(model.group,trueX)
tab.group[isimu,1] <- quality$TTPR
tab.group[isimu,2] <- quality$TFPR
tab.group[isimu,3] <- quality$TD
  
  
#res.group <- select.sgspls(model.group)
#res.group$group.size.X
#res.group$group.size.Y
  
#### Sparse Group PLS
res.tun.sgpls <- choice.tuning.sgpls.X.comp(X,Y,folds = folds,ncomp=1,grid.gX,grid.alpha.X,setseed=choicesetseed,ind.block.x=ind.block.x,upper.lambda=1000000000)
model.sparse <- Sparse.Group.spls.BP(X,Y,ncomp=1,mode="regression",max.iter=500,tol=1e-06,keepX=rep(res.tun.sgpls$keepX,1),ind.block.x=ind.block.x,alpha.x=res.tun.sgpls$alphaX,upper.lambda=1000000000)

if(norm(model.sparse$loadings$X[,1]-theta.xx1) < norm(-model.sparse$loadings$X[,1]-theta.xx1)) temp <- model.sparse$loadings$X[,1] else  temp <- -model.sparse$loadings$X[,1]
res.sparseX1 <- cbind(res.sparseX1,temp)

#res.sparse <- select.sgspls(model.sparse)
#res.sparse$group.size.X
#res.sparse$group.size.Y

quality <- quality.simu.X(model.sparse,trueX)
tab.sparse[isimu,1] <- quality$TTPR
tab.sparse[isimu,2] <- quality$TFPR
tab.sparse[isimu,3] <- quality$TD


#### SPLS

#res.tun.spls <- choice.tuning.spls.X.comp(X,Y,folds = folds,ncomp=1,grid.X,setseed=choicesetseed)
#res.tun.spls.mix <- choice.tuning.spls.mix.X.comp(X,Y,folds = folds,ncomp=1,grid.X,setseed=choicesetseed)

model.spls <- spls.BP(X,Y,ncomp=1,mode="regression",max.iter=500,tol=1e-06,keepX=rep(60,1))

if(norm(model.spls$loadings$X[,1]-theta.xx1) < norm(-model.spls$loadings$X[,1]-theta.xx1)) temp <- model.spls$loadings$X[,1] else  temp <- -model.spls$loadings$X[,1]
res.splsX1 <- cbind(res.splsX1,temp)

quality <- quality.simu.X(model.spls,trueX)
tab.spls[isimu,1] <- quality$TTPR
tab.spls[isimu,2] <- quality$TFPR
tab.spls[isimu,3] <- quality$TD

#res.spls <- select.spls(model.spls)
#res.spls$select.X
#res.spls$select.Y
}

res.groupXm1 <- rowMeans(res.groupX1)
res.splsXm1 <- rowMeans(res.splsX1)
res.sparseXm1 <- rowMeans(res.sparseX1)




res.qua.group <- rbind(res.qua.group,colMeans(tab.group))
res.qua.sparse <- rbind(res.qua.sparse,colMeans(tab.sparse))
res.qua.spls <- rbind(res.qua.spls,colMeans(tab.spls))


res.tab <- rbind(colMeans(tab.spls),colMeans(tab.group),colMeans(tab.sparse))
colnames(res.tab) <- c("TTPR","TFPR","TD")
rownames(res.tab) <- c("sPLS","gPLS","sgPLS")
res.tab <- signif(res.tab,3)

name.tab <-paste("tab-simu1-",kk,".tex",sep="")
print(xtable(res.tab,digits=4),file=paste("Dropbox/PROJECT-A-STAR-PLS/SIMULATION/SIMU-KEEP-corr/",name.tab,sep=""))

x <- 1:400
pdf(file=paste("Dropbox/PROJECT-A-STAR-PLS/SIMULATION/SIMU-KEEP-corr/simu1u1-",kk,".pdf",sep=""),width=8,height=11)
par(mfrow=c(2,2))
plot(theta.xx1~x,ylab=c(""),xlab="",col="blue")
legend("topleft",legend="Original loading: u1",col="blue",pch=1)
plot(res.groupXm1~x,col="red",ylab=c(""),xlab="")
legend("bottomleft",legend="Group PLS loading: u1",col="red",pch=1)
plot(res.splsXm1~x,col="red",ylab=c(""),xlab="")
legend("bottomleft",legend="sPLS loading: u1",col="red",pch=1)
plot(res.sparseXm1~x,col="red",ylab=c(""),xlab="")
legend("bottomleft",legend="Group sPLS loading: u1",col="red",pch=1)
dev.off()

}
save.image(file="Dropbox/PROJECT-A-STAR-PLS/SIMULATION/SIMU-KEEP-corr/simulation1.RData")


if(FALSE){
pdf(file="Dropbox/PROJECT-A-STAR-PLS/SIMULATION/simu7-noise.pdf",width=8,height=11)
par(mfrow=c(2,1))
plot(res.qua.group[,1]~set.sigma.e,type="b",ylim=c(0.5,1),col="blue",ylab="TTPR",xlab="Noise (standard deviation)")
lines(res.qua.spls[,1]~set.sigma.e,type="b",ylim=c(0.5,1),col="red")
lines(res.qua.sparse[,1]~set.sigma.e,type="b",ylim=c(0.5,1),col="green")
legend("bottomleft",legend=c("sPLS","gPLS","sgPLS"),lty=1,col=c("red","blue","green"),pch=1)

plot(res.qua.group[,3]~set.sigma.e,ylim=c(0,max(res.qua.group[,3],res.qua.spls[,3],res.qua.sparse[,3])),type="b",col="blue",ylab="TD",xlab="Noise (standard deviation)")
lines(res.qua.spls[,3]~set.sigma.e,type="b",col="red")
lines(res.qua.sparse[,3]~set.sigma.e,type="b",col="green")
legend("bottomright",legend=c("sPLS","gPLS","sgPLS"),lty=1,col=c("red","blue","green"),pch=1)

dev.off()


x <- 1:400
pdf(file="Dropbox/PROJECT-A-STAR-PLS/SIMULATION/simu7u1.pdf",width=8,height=11)
par(mfrow=c(2,2))
plot(theta.xx1~x,ylab=c(""),xlab="",col="blue")
legend("topleft",legend="Original loading: u1",col="blue",pch=1)
plot(res.groupXm1~x,col="red",ylab=c(""),xlab="")
legend("bottomleft",legend="Group PLS loading: u1",col="red",pch=1)
plot(res.splsXm1~x,col="red",ylab=c(""),xlab="")
legend("bottomleft",legend="sPLS loading: u1",col="red",pch=1)
plot(res.sparseXm1~x,col="red",ylab=c(""),xlab="")
legend("bottomleft",legend="Group sPLS loading: u1",col="red",pch=1)
dev.off()

}


