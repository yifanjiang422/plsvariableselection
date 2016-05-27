
data.construct <- function(n,p,q,off.var,ncomp,plevel,qlevel,rho){
  data.raw.c.list <- list()
  data.raw.com <- matrix(rnorm(n*ncomp,0,1),n,ncomp)          #generate the matrix with random variables from normdis with size n-off.var by n
  
  data.raw.c <- matrix(c(sample(c(-1,1,1.5),p,replace=TRUE))
                       ,ncomp,p,byrow=TRUE)                   #
  index <- seq(1:p)
  zero.index <- sample(3:p,off.var)
  data.raw.c[1,zero.index]=0
  if(ncomp > 1){
    zero.index2 <- index[-zero.index]
    
    data.raw.c[2,zero.index2]=0
  }
  
  data.raw.d <- matrix(rep(runif(q,0.5,2),ncomp),ncomp,q)
  d.index <- seq(1:q)
  sign.d <- rep(c(1,-1),(length(d.index)/2))
  if(ncomp>1){
    data.raw.d[2,d.index]<- sign.d*data.raw.d[1,rev(d.index)]
  }
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
  
  data.raw.x <- data.raw.com %*% data.raw.c + (data.raw.errx)
  data.raw.z <- data.raw.com %*% data.raw.d + (data.raw.erry)
  
  x.cat <- catConvert(data.raw.x,plevel)
  z.cat <- catConvert(data.raw.z,qlevel)
  
  if(ncomp>1){
    zero.index.sum=c(zero.index,zero.index2)
  }else{
    zero.index.sum=zero.index
  }
  
  return(list(x.cat=x.cat,z.cat=z.cat,
              data.raw.x=data.raw.x, data.raw.z=data.raw.z,
              zero.index=zero.index.sum,
              data.raw.c=data.raw.c,
              data.raw.d=data.raw.d,
              data.raw.com=data.raw.com))
}