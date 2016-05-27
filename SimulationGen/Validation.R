valid.count <- function(sim.no, mode, loadings, off.var, ncomp, zero.index, plevel,qlevel){
  invalid.count = c()                                                  #Container for inaccurate capture of insignificant variables  
  m.invalid = c()                                                      #Container to report final num of invalid capture of variables
  acc.tpr=c()
  acc.tpr.temp=c()
  acc.td=c()
  acc.td.temp=c()
  test=c()
  m.spars.calc=list()
  m1.index=NULL
  if (mode == "PLSCA"){
    spars.group.ind.ex=c()                                             #spars group member index defined for the simulation. Known
    spars.group.ind <- (zero.index * plevel[1]-plevel[1]+1)
    for(i in 1:length(spars.group.ind)){
      spars.group.ind.ex <- c(spars.group.ind.ex,seq(from=(spars.group.ind[i]),to=(spars.group.ind[i]+plevel[1]-1),by=1))
    }
    calc.spars.ind=which(loadings$X == 0)%%(plevel[1] * length(plevel))  #Spars group calculated from simulation
    calc.spars.ind[which(calc.spars.ind==0)]=(plevel[1] * length(plevel))#revert 0 to off.var*length(plevel)
    if(length(calc.spars.ind) <ncomp* off.var * plevel[1]){
      calc.spars.ind = c(calc.spars.ind,rep(0,((ncomp*off.var * plevel[1])-length(calc.spars.ind))))
    }
    for(i in 1:(length(calc.spars.ind)-1)){
      if(calc.spars.ind[i]>calc.spars.ind[i+1]){
        m1.index=i
        break
      }
    }
    if(is.null(m1.index)){
      if(ncomp==1){
        m.spars.calc[[1]]=calc.spars.ind
      }
      if(ncomp>1){
        m.spars.calc[[1]]=calc.spars.ind[1:length(which(loadings$X[,1] == 0))]
        m.spars.calc[[2]]=calc.spars.ind[(length(which(loadings$X[,1] == 0))+1):length(calc.spars.ind)]
      }
    }else{
      m.spars.calc[[1]]=calc.spars.ind[1:m1.index]
      m.spars.calc[[2]]=calc.spars.ind[(m1.index+1):length(calc.spars.ind)]
    }
    
    #if more zero index exists in calculated result, 
    if(length(calc.spars.ind)>(length(zero.index)*plevel[1])){
      fpno=(length(calc.spars.ind)-length(zero.index)*plevel[1])/plevel[1]
    }else{
      fpno=0
    }
    
    m.spars.ex <- matrix(spars.group.ind.ex,nrow= ncomp,
                         ncol=length(spars.group.ind.ex)/ncomp,
                         byrow=TRUE)
    m.full.ind = seq(from=1,to=plevel[1]*length(plevel))
    
    for(i in 1:ncomp){
      m.spars.ex[i,]=sort(as.numeric(m.spars.ex[i,]))
      
      for(j in 1:ncomp){
        temp = off.var*plevel[1]-length(intersect(m.spars.ex[i,],m.spars.calc[[j]]))        #comparing the location of 0s in each column of loadings and record the mismatches
        invalid.count = c(invalid.count,temp)
        
        m.tp.ind = setdiff(m.full.ind,m.spars.ex[i,])
        m.tp.ind.calc = setdiff(m.full.ind,m.spars.calc[[j]])
        m.fn.ind.calc = setdiff(m.spars.calc[[j]],m.spars.ex[i,])
        acc.tpr.temp=c(acc.tpr.temp,length(intersect(m.tp.ind,m.tp.ind.calc))/(length(intersect(m.tp.ind,m.tp.ind.calc))+length(m.fn.ind.calc)))
        
        m.fp.ind.calc = setdiff(m.tp.ind.calc,m.tp.ind)
        acc.td.temp=c(acc.td.temp,length(m.fp.ind.calc)+length(m.fn.ind.calc))
 
      }
      m.invalid=c(m.invalid,min(invalid.count))
      invalid.count=c()
      
      acc.tpr=c(acc.tpr,max(acc.tpr.temp))
      acc.tpr.temp=c()
      
      acc.td=c(acc.td,min(acc.td.temp))
      acc.td.temp=c()
      
    }
    tpr=mean(acc.tpr)
    td=sum(acc.td)
    acc.per <- 1-((sum(m.invalid)/plevel[1])+fpno)/(off.var*ncomp)                     #Calculate the accurate percentage of the program
  }
  if (mode == "PLSQT"){
    m.zero.index <- matrix(zero.index,ncomp,off.var,byrow=TRUE)
    m.zero.index <- t(apply(m.zero.index,1,sort))
    spars.calc.ind <- which(loadings$X==0)%%(length(plevel))
    spars.calc.ind[which(spars.calc.ind==0)]=length(plevel)
    m.full.ind = seq(from=1,to=length(plevel))
    
    for(i in 1:(length(spars.calc.ind)-1)){
      if(spars.calc.ind[i]>spars.calc.ind[i+1]){
        m1.index=i
        break
      }
    }
    if(is.null(m1.index)){
      if(ncomp==1){
        m.spars.calc[[1]]=spars.calc.ind
      }
      if(ncomp>1){
        m.spars.calc[[1]]=spars.calc.ind[1:length(which(loadings$X[,1] == 0))]
        m.spars.calc[[2]]=spars.calc.ind[(length(which(loadings$X[,1] == 0))+1):length(spars.calc.ind)]
      }
    }else{
      m.spars.calc[[1]]=spars.calc.ind[1:m1.index]
      m.spars.calc[[2]]=spars.calc.ind[(m1.index+1):length(spars.calc.ind)]
    }
    
    if(length(spars.calc.ind)>length(zero.index)){
      fpno=length(spars.calc.ind)-length(zero.index)
    }else{
      fpno=0
    }
    
    #m.spars.calc <- matrix(spars.calc.ind,ncomp,m1.index,byrow=TRUE)
    for(i in 1:ncomp){
      for(j in 1:ncomp){
        temp = off.var - length(intersect(m.spars.calc[[i]],m.zero.index[j,]))
        invalid.count = c(invalid.count,temp)
        
        m.tp.ind = setdiff(m.full.ind,m.zero.index[i,])
        m.tp.ind.calc = setdiff(m.full.ind,m.spars.calc[[j]])
        m.fn.ind.calc = setdiff(m.spars.calc[[j]],m.zero.index[i,])
        acc.tpr.temp=c(acc.tpr.temp,length(intersect(m.tp.ind,m.tp.ind.calc))/(length(intersect(m.tp.ind,m.tp.ind.calc))+length(m.fn.ind.calc)))
        
        m.fp.ind.calc = setdiff(m.tp.ind.calc,m.tp.ind)
        acc.td.temp=c(acc.td.temp,length(m.fp.ind.calc)+length(m.fn.ind.calc))
        
        
      }
      m.invalid=c(m.invalid,min(invalid.count))
      invalid.count=c()
      
      acc.tpr=c(acc.tpr,max(acc.tpr.temp))
      acc.tpr.temp=c()
      
      acc.td=c(acc.td,min(acc.td.temp))
      acc.td.temp=c()
    }
    acc.per <- 1-((sum(m.invalid)+fpno)/(off.var*ncomp))
    tpr=mean(acc.tpr)
    td=sum(acc.td)
  }
  return(list(invsim.no=test,invalid.v=sum(m.invalid),accuracy=acc.per,td=td,tpr=tpr))
}