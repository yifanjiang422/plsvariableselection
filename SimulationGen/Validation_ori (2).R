valid.count <- function(sim.no, mode, loadings, off.var, ncomp, zero.index, plevel,qlevel){
  invalid.count = c()                                                  #Container for inaccurate capture of insignificant variables  
  m.invalid = c()                                                      #Container to report final num of invalid capture of variables
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
      calc.spars.ind = c(calc.spars.ind,rep(100,((ncomp*off.var * plevel[1])-length(calc.spars.ind))))
    }
    for(i in 1:(length(calc.spars.ind)-1)){
      if(calc.spars.ind[i]>calc.spars.ind[i+1]){
        m1.index=i
        break
      }
    }
    if(is.null(m1.index)){
      m.spars.calc[[1]]=calc.spars.ind
    }else{
      m.spars.calc[[1]]=calc.spars.ind[1:m1.index]
      m.spars.calc[[2]]=calc.spars.ind[(m1.index+1):length(calc.spars.ind)]
    }
    
    if(length(calc.spars.ind)>(length(zero.index)*plevel[1])){
      fpno=(length(calc.spars.ind)-length(zero.index)*plevel[1])/plevel[1]
    }else{
      fpno=0
    }
    
    m.spars.ex <- matrix(spars.group.ind.ex,nrow= ncomp,
                         ncol=length(spars.group.ind.ex)/ncomp,
                         byrow=TRUE)
#     if(length(calc.spars.ind)>(ncomp*off.var * plevel[1])){
#       fpno=length(calc.spars.ind)-(ncomp*off.var * plevel[1])
#       calc.spars.ind=calc.spars.ind[1:(ncomp*off.var * plevel[1])]
#     }else{
#       fpno=0
#     }
#     for(vindex in 1:(ncomp*off.var * plevel[1]-1)){
#       if(calc.spars.ind[vindex]>calc.spars.ind[vindex+1]){
#         m1.index=vindex
#         break
#       }else{
#         m1.index=(length(calc.spars.ind)/ncomp)
#       }
#     }
#    m.spars.calc <- matrix(calc.spars.ind,ncomp,
#                           m1.index,
#                           byrow=TRUE)
    
    for(i in 1:ncomp){
      m.spars.ex[i,]=sort(as.numeric(m.spars.ex[i,]))
      for(j in 1:ncomp){
        temp = off.var*plevel[1]-length(intersect(m.spars.ex[i,],m.spars.calc[[j]]))        #comparing the location of 0s in each column of loadings and record the mismatches
        invalid.count = c(invalid.count,temp)
      }
      m.invalid=c(m.invalid,min(invalid.count))
      invalid.count=c()
    }
    
    acc.FPR <- sum(m.invalid)
    acc.per <- 1-((sum(m.invalid)/plevel[1])+fpno)/(off.var*ncomp)                     #Calculate the accurate percentage of the program
  }
  if (mode == "PLSQT"){
    m.zero.index <- matrix(zero.index,ncomp,off.var,byrow=TRUE)
    m.zero.index <- t(apply(m.zero.index,1,sort))
    spars.calc.ind <- which(res.qt.vs$loadings$X==0)%%(length(plevel))
    spars.calc.ind[which(spars.calc.ind==0)]=length(plevel)
    
    
    for(i in 1:(length(spars.calc.ind)-1)){
      if(spars.calc.ind[i]>spars.calc.ind[i+1]){
        m1.index=i
        break
      }
    }
    if(is.null(m1.index)){
      m.spars.calc[[1]]=spars.calc.ind
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
      }
      m.invalid=c(m.invalid,min(invalid.count))
      invalid.count=c()
    }
    acc.per <- 1-((sum(m.invalid)+fpno)/(off.var*ncomp))
  }
  return(list(invsim.no=test,invalid.v=sum(m.invalid),accuracy=acc.per))
}