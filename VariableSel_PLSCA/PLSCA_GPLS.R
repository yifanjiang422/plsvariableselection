plsca_gpls <- function(n,p,q,off.var,ncomp,plevel,qlevel,rho,data1,data2){
  
  #create ind.block.x and ind.block.z
  ind.block.x<-ind.block.z<-NULL
  ind.block.x<-plevel[1]
  ind.block.z<-qlevel[1]
  for(i in 1:(p-2)){
    ind.block.x=cbind(ind.block.x,sum(plevel[1:(i+1)]))
  }
  for(i in 1:(q-2)){
    ind.block.z=cbind(ind.block.z,sum(qlevel[1:(i+1)]))
  }  #end create ind.block.x and z
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
  res.ca.vs <- plscaGPLS(data1, data2, ncomp = ncomp, mode = "regression", keepX = (rep(p-off.var,p)), 
                         keepZ = (rep(q,q)), ind.block.x = ind.block.x
                         ,ind.block.z = ind.block.z, masses = masses, weights=weights)
  return(res.ca.vs)
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
}