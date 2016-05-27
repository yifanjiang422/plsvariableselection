catConvert <- function(raw,level){
  var.no <- ncol(raw)
  var.ob <- nrow(raw)
  temp <- vector()
  #raw<-round(raw,4)
  raw<-apply(raw,2,as.character)
  for(i in 1:var.no){
   lval <- (max(as.numeric(raw[,i]))-min(as.numeric(raw[,i])))/level[i]
   temp <- raw[,i]
   temp.num <- as.numeric(raw[,i])
   for(j in 1:level[i]){
     if(j==level[i]){
       ind=j
       j=j+j
     }
     else{
       ind=j
     }
     temp[which(temp.num>=(max(as.numeric(raw[,i]))-j*lval))]<-paste(i,".",ind,sep="")
     temp.num[which(temp.num>=(max(as.numeric(raw[,i]))-j*lval))]<-min(as.numeric(raw[,i]))-(level[i]+1)*lval
   }  
   raw[,i]<-temp
  }
  return(raw)
  
  
}
#test1[which(as.numeric(test1)>5)]<-paste("V", i, "L", i)