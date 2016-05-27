dummy_GPLS <- function(n,p,q,off.var,ncomp,plevel,qlevel,rho,data1,data2){

  #create ind.block.x and ind.block.z
  ind.block.x<-ind.block.z<-NULL
  ind.block.x<-plevel[1]
  ind.block.z<-qlevel[1]
  for(i in 1:(p-2)){
    ind.block.x=cbind(ind.block.x,sum(plevel[1:(i+1)]))
  }
  for(i in 1:(q-2)){
    ind.block.z=cbind(ind.block.z,sum(qlevel[1:(i+1)]))
  }
  #data1[,(ind.block.x[p-2]+1):dim(data1)[2]]=data1[,(ind.block.x[p-2]+1):dim(data1)[2]]
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
  res.dum.vs <- dummyGPLS(data1, data2, ncomp = ncomp, mode = "regression", keepX = (rep(p-off.var,p)), 
                        keepZ = (rep(q,q)), ind.block.x = ind.block.x
                        ,ind.block.z = ind.block.z, masses = masses, weights=weights)

  return(res.dum.vs)
}