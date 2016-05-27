quant_spls <- function(n,p,q,off.var,ncomp,plevel,qlevel,rho,data1,data2){

  #create ind.block.x and ind.block.z
#   ind.block.x<-ind.block.z<-NULL
#   ind.block.x<-plevel[1]
#   ind.block.z<-qlevel[1]
#   for(i in 1:(p-off.var)){
#     ind.block.x=cbind(ind.block.x,sum(plevel[1:(i+1)]))
#   }
#   for(i in 1:(q-2)){
#     ind.block.z=cbind(ind.block.z,sum(qlevel[1:(i+1)]))
#   }
  
  #end create ind.block.x and z
  
  quantscale=sqrt((n-1)/n) #correction factor for optimum scaling method
  Quant.data=osquant(data1,data2,ncomp,quantscale)
  
  res.qt.vs <- qspls(Quant.data$QuantX, Quant.data$QuantZ, ncomp=ncomp, keepX = (rep(p-off.var,ncomp)), keepY = (rep(q,ncomp)),
                     U=Quant.data$U, V=Quant.data$V, Omega=Quant.data$Omega, Xi=Quant.data$Xi, noiteration=Quant.data$noiteration)
  return(res.qt.vs)
}

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
