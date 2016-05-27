valid.count <- function(sim.no, mode, loadings, off.var, ncomp, zero.index, plevel,qlevel){
  invalid.count = c()                                                  #Container for inaccurate capture of insignificant variables  
  m.invalid = c()                                                      #Container to report final num of invalid capture of variables
  test=c()
  if (mode == "PLSCA"){
    spars.group.ind.ex=c()                                             #spars group member index defined for the simulation. Known
    spars.group.ind <- (zero.index * plevel[1]-plevel[1]+1)
    for(i in 1:length(spars.group.ind)){
      spars.group.ind.ex <- c(spars.group.ind.ex,seq(from=(spars.group.ind[i]),to=(spars.group.ind[i]+plevel[1]-1),by=1))
    }
    calc.spars.ind=which(loadings$X == 0)%%(off.var * length(plevel))  #Spars group calculated from simulation
    calc.spars.ind[which(calc.spars.ind==0)]=(off.var * length(plevel))#revert 0 to off.var*length(plevel)
    m.spars.ex <- matrix(spars.group.ind.ex,nrow= ncomp,
                         ncol=length(spars.group.ind.ex)/ncomp,
                         byrow=TRUE)
    m.spars.calc <- matrix(calc.spars.ind,ncomp,
                           (length(calc.spars.ind)/ncomp),
                           byrow=TRUE)
    
    for(i in 1:ncomp){
      m.spars.ex[i,]=sort(as.numeric(m.spars.ex[i,]))
      for(j in 1:ncomp){
        temp = length(which(m.spars.ex[j,]!=m.spars.calc[j,]))        #comparing the location of 0s in each column of loadings and record the mismatches
        invalid.count = c(invalid.count,temp)
        if(length(m.spars.ex[j,])!=length(m.spars.calc[j,])){
          test =c(test,sim.no)
        }
      }
      m.invalid=c(m.invalid,min(invalid.count))
      invalid.count=c()
    }
    acc.per <- 1-(sum(m.invalid)/plevel[1])/(off.var*ncomp)                     #Calculate the accurate percentage of the program
  }
  if (mode == "PLSQT"){
    m.zero.index <- matrix(zero.index,ncomp,off.var,byrow=TRUE)
    m.zero.index <- t(apply(m.zero.index,1,sort))
    spars.calc.ind <- which(res.qt.vs$loadings$X==0)%%(length(plevel))
    spars.calc.ind[which(spars.calc.ind==0)]=length(plevel)
    m.spars.calc <- matrix(spars.calc.ind,ncomp,off.var,byrow=TRUE)
    for(i in 1:ncomp){
      for(j in 1:ncomp){
        temp = length(which(m.spars.calc[i,] != m.zero.index[j,]))
        invalid.count = c(invalid.count,temp)
      }
      m.invalid=c(m.invalid,min(invalid.count))
      invalid.count=c()
    }
    acc.per <- 1-(sum(m.invalid)/(off.var*ncomp))
  }
  return(list(invsim.no=test,invalid.v=m.invalid,accuracy=acc.per))
}