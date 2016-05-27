
###########Loading Resources############
library(MASS)
library(TExPosition)
library(ExPosition)
library(ggplot2)
library(reshape2)
library(spls)
library(binr)
library(sgPLS)
library(mvtnorm)
library(Matrix)
source("SimulationGen/catConvert.R")
source("SimulationGen/dataconstruct.R")
source("SimulationGen/Validation.R")
source("PLSR/GPLSos.R")
source("CatPLS/createIndMatrix.R")
source("VariableSel_Quant/InternalFunc/Quantification.R")
source("VariableSel_Quant/InternalFunc/Quant_spls.R")
source("VariableSel_Quant/InternalFunc/Quant_spls_step2.R")
source("VariableSel_Quant/Quant_SPLS.R")
source("VariableSel_PLSCA/PLSCA_GPLS.R")
source("VariableSel_PLSCA/InternalFunc/plscaGPLS.R")
source("VariableSel_PLSCA/InternalFunc/plscaGPLSStep1.R")
source("VariableSel_PLSCA/InternalFunc/plscaGPLSStep2.R")
source("VariableSel_PLSCA/InternalFunc/rowNormsvs.R")
source("VariableSel_PLSCA/InternalFunc/caSVD.R")
source("VariableSel_PLSCA/InternalFunc/makeRowProfilesvs.R")
source("VariableSel_PLSCA/InternalFunc/caNormvs.R")
source("VariableSelDummy/InternalFunc/dummyGPLS.R")
source("VariableSelDummy/InternalFunc/dummyGPLSStep1.R")
source("VariableSelDummy/dummy_GPLS.R")
###########End Loading Resources############

###########Define User Input############
res.qt.res=c()
res.dum.res=c()
res.ca.res=c()
n <- 20                                        #Define size of observation
p <- 10                                        #No. of X variables
q <- 4                                          #No. of Z variables (need to be even number at the moment)
v.offvar <- c(6,7,8)
for(off.var in v.offvar){
#off.var <- 60                                   #No. of variables not relevant
ncomp <- 1                                      #No. of Components to be seletected
set.seed(13)                                    #determine the levels in each cat variable randomly (change in dataconstruct as well)
plevel <- rep(2,p)                              #Levels in X 
qlevel <- rep(2,q)                              #Levels in Z
#plevel <- sample(seq(3,4,by=1),p,replace=TRUE) #Generate different no. of levels for X variable
#qlevel <- sample(seq(3,4,by=1),q,replace=TRUE) #Generate different no. of levels for X variable
rho=0.5
simulation.no=100
###########End User Input############
res.ca.acc=c()
res.ca.invalid.ind=matrix(,simulation.no,ncomp)
res.ca.invsim.no=c()
res.dum.acc=c()
res.dum.invalid.ind=matrix(,simulation.no,ncomp)
res.dum.invsim.no=c()
res.qt.acc=c()
res.qt.invalid.ind=matrix(,simulation.no,ncomp)

###########Built Data Matrix##############

for(i in 1:simulation.no){
  if(i==38){
    i=38
  }
  data.list<-data.construct(n,p,q,off.var,ncomp,plevel,qlevel,rho)
  
  ###########End Build######################
  
  ###########PLSCA########
  plsca.data1=createIndMatrix(data.list$x.cat,off.var,plevel)
  plsca.data2=createIndMatrix(data.list$z.cat,0,qlevel)
  res.ca.vs <- plsca_gpls(n,p,q,off.var,ncomp,plevel,qlevel,rho,plsca.data1,plsca.data2)
  ###########PLSCAEND######
  
  ##########QuantPLS#######
  quant.data1=as.matrix(data.list$x.cat)
  quant.data2=as.matrix(data.list$z.cat)
  res.qt.vs <- quant_spls(n,p,q,off.var,ncomp,plevel,qlevel,rho,quant.data1,quant.data2)
  ##########QuantPLSEND####
  
  ##########Dummy##########
  dummy.data1=createIndMatrix(data.list$x.cat,off.var,plevel)
  dummy.data2=createIndMatrix(data.list$z.cat,0,qlevel)
  res.dum.vs <- dummy_GPLS(n,p,q,off.var,ncomp,plevel,qlevel,rho,dummy.data1,dummy.data2)
  
  ###############Data Validation section################
  res.ca.valid <- valid.count(sim.no = i, mode = "PLSCA", loadings = res.ca.vs$loadings, off.var=off.var, ncomp=ncomp, zero.index = data.list$zero.index, 
                              plevel=plevel, qlevel=qlevel)
  res.ca.acc=c(res.ca.acc,res.ca.valid$accuracy)
  res.ca.invalid.ind[i,]=res.ca.valid$invalid.v
  if(!is.null(res.ca.valid$invsim.no)){
    res.ca.invsim.no = c(res.ca.invsim.no,res.ca.valid$invsim.no)
  }
  
  res.dum.valid <- valid.count(sim.no = i, mode = "PLSCA", loadings = res.dum.vs$loadings, off.var=off.var, ncomp=ncomp, zero.index = data.list$zero.index, 
                               plevel=plevel, qlevel=qlevel)
  res.dum.acc=c(res.dum.acc,res.dum.valid$accuracy)
  res.dum.invalid.ind[i,]=res.dum.valid$invalid.v
  if(!is.null(res.dum.valid$invsim.no)){
    res.dum.invsim.no = c(res.dum.invsim.no,res.dum.valid$invsim.no)
  }
  
  res.qt.valid <- valid.count(sim.no = i, mode = "PLSQT", loadings = res.qt.vs$loadings, off.var=off.var, ncomp=ncomp, zero.index = data.list$zero.index, 
                              plevel=plevel, qlevel=qlevel)
  res.qt.acc=c(res.qt.acc,res.qt.valid$accuracy)
  res.qt.invalid.ind[i,]=res.qt.valid$invalid.v
  
  print(i)
}
res.qt.res <- c(res.qt.res,max(res.qt.acc))
res.qt.res <- c(res.qt.res,mean(res.qt.acc))
res.qt.res <- c(res.qt.res,min(res.qt.acc))

res.ca.res <- c(res.ca.res,max(res.ca.acc))
res.ca.res <- c(res.ca.res,mean(res.ca.acc))
res.ca.res <- c(res.ca.res,min(res.ca.acc))

res.dum.res<- c(res.dum.res,max(res.dum.acc))
res.dum.res<- c(res.dum.res,mean(res.dum.acc))
res.dum.res<- c(res.dum.res,min(res.dum.acc))
}
