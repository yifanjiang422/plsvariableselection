catConvert <- function(rawx,level){
  var.no <- ncol(rawx)
  var.ob <- nrow(rawx)
  rownames(rawx)<-paste(1:var.ob,sep="")
  for(i in 1:var.no){
   sec.bin <- bins.quantiles(rawx[,i],level[i],level[i])
   sec.lo.ind <- sec.bin$binlo
   sec.hi.ind <- sec.bin$binhi
   temp.v<-as.matrix(sort(rawx[,i]))
   for(j in 1:level[i]){
     temp.v[sec.lo.ind[j]:sec.hi.ind[j],1]<-as.numeric(i+(j/100))
   }  
   rawx[,i]<-as.matrix(temp.v[order(rownames(temp.v))])
  }
  return(rawx)
  
  
}
#test1[which(as.numeric(test1)>5)]<-paste("V", i, "L", i)