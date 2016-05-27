
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
source("VariableSel_Quant/InternalFunc/Quant_spls_step1.R")
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
source("MixdataType/catConvert.R")
source("Mixdatatype/GPLSos.R")
source("Mixdatatype/createIndMatrix.R")
source("MixdataType/InternalFunc_Quant/Quant_spls.R")
source("MixdataType/InternalFunc_Quant/Quant_spls_step2.R")
source("MixdataType/mixQuant_SPLS.R")

###########End Loading Resources############

###########Define User Input############
res.qt.res=c()
res.dum.res=c()
res.ca.res=c()
res.mix.res=c()
res.mean.tpr=c()
res.sum.td=c()
n <- 100                                       #Define size of observation
p <- 100                                        #No. of X variables
q <- 4                                          #No. of Z variables (need to be even number at the moment)
#v.offvar <- c(45,30,20,10,5)
v.offvar <- c(90,70,50,30,10)
for(off.var in v.offvar){
#off.var <- 60                                   #No. of variables not relevant
ncomp <- 2                                      #No. of Components to be seletected
set.seed(25)                                    #determine the levels in each cat variable randomly (change in dataconstruct as well)
plevel <- rep(3,p)                              #Levels in X 
qlevel <- rep(3,q)                              #Levels in Z
#plevel <- sample(seq(3,4,by=1),p,replace=TRUE) #Generate different no. of levels for X variable
#qlevel <- sample(seq(3,4,by=1),q,replace=TRUE) #Generate different no. of levels for X variable
rho=0.5
simulation.no=100
###########End User Input############
res.ca.acc=c()
res.ca.invalid.ind=matrix(,simulation.no,ncomp)
res.ca.invsim.no=c()
res.ca.tpr=c()
res.ca.td=c()
res.dum.acc=c()
res.dum.invalid.ind=matrix(,simulation.no,ncomp)
res.dum.invsim.no=c()
res.dum.tpr=c()
res.dum.td=c()
res.qt.acc=c()
res.qt.invalid.ind=matrix(,simulation.no,ncomp)
res.qt.tpr=c()
res.qt.td=c()
res.mix.acc=c()
res.mix.invalid.ind=matrix(,simulation.no,ncomp)
res.mix.tpr=c()
res.mix.td=c()
res.spls.tpr=c()
res.spls.td=c()

###########Built Data Matrix##############

for(i in 1:simulation.no){
  if(i==17){
    i=17
  }
  data.list<-data.construct(n,p,q,off.var,ncomp,plevel,qlevel,rho)
  #mixdata.list<-mixdata.construct(n,p,q,off.var,ncomp,plevel,qlevel,rho)
  ###########End Build######################
  
  ###########SPLS#########
  spls.data1<-data.list$data.raw.x
  spls.data2<-data.list$data.raw.z
  res.spls.vs<-sPLS(X=spls.data1,Y=spls.data2,ncomp=ncomp,mode="regression",keepX=rep(p-off.var,2))
  
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
  
  ########MIX#############
  quant.data1.cat=as.matrix(data.list$mixx.cat)
  quant.data2.cat=as.matrix(data.list$mixz.cat)
  quant.data1.qnt=as.matrix(data.list$mixx.quant)
  quant.data2.qnt=as.matrix(data.list$mixz.quant)
  res.mix.vs <- mixquant_spls(n,p,q,off.var,ncomp,plevel,qlevel,rho,quant.data1.cat,quant.data2.cat,quant.data1.qnt,quant.data2.qnt)
  #######MixEnd##########
  
  
  res.ca.valid <- valid.count(sim.no = i, mode = "PLSCA", loadings = res.ca.vs$loadings, off.var=off.var, ncomp=ncomp, zero.index = data.list$zero.index, 
                             plevel=plevel, qlevel=qlevel)
  res.ca.acc=c(res.ca.acc,res.ca.valid$accuracy)
  res.ca.invalid.ind[i,]=res.ca.valid$invalid.v
  res.ca.td=c(res.ca.td,res.ca.valid$td)
  res.ca.tpr=c(res.ca.tpr,res.ca.valid$tpr)
  ##############################################################
  res.dum.valid <- valid.count(sim.no = i, mode = "PLSCA", loadings = res.dum.vs$loadings, off.var=off.var, ncomp=ncomp, zero.index = data.list$zero.index, 
                               plevel=plevel, qlevel=qlevel)
  res.dum.acc=c(res.dum.acc,res.dum.valid$accuracy)
  res.dum.invalid.ind[i,]=res.dum.valid$invalid.v
  res.dum.td=c(res.dum.td,res.dum.valid$td)
  res.dum.tpr=c(res.dum.tpr,res.dum.valid$tpr)
  ###############################################################
  res.qt.valid <- valid.count(sim.no = i, mode = "PLSQT", loadings = res.qt.vs$loadings, off.var=off.var, ncomp=ncomp, zero.index = data.list$zero.index, 
                              plevel=plevel, qlevel=qlevel)
  res.qt.acc=c(res.qt.acc,res.qt.valid$accuracy)
  res.qt.invalid.ind[i,]=res.qt.valid$invalid.v
  res.qt.td=c(res.qt.td,res.qt.valid$td)
  res.qt.tpr=c(res.qt.tpr,res.qt.valid$tpr)
  ##############################################################
  res.mix.valid <- valid.count(sim.no = i, mode = "PLSQT", loadings = res.mix.vs$loadings, off.var=off.var, ncomp=ncomp, zero.index = data.list$zero.index, 
                                  plevel=plevel, qlevel=qlevel)
  res.mix.acc=c(res.mix.acc,res.mix.valid$accuracy)
  res.mix.invalid.ind[i,]=res.mix.valid$invalid.v
  res.mix.td=c(res.mix.td,res.mix.valid$td)
  res.mix.tpr=c(res.mix.tpr,res.mix.valid$tpr)
  #############################################################
  res.spls.valid<-valid.count(sim.no = i, mode = "PLSQT", loadings = res.spls.vs$loadings, off.var=off.var, ncomp=ncomp, zero.index = data.list$zero.index, 
                              plevel=plevel, qlevel=qlevel)
  res.spls.td=c(res.spls.td,res.spls.valid$td)
  res.spls.tpr=c(res.spls.tpr,res.spls.valid$tpr)
  
  
  print(i)
}
res.mean.tpr<-c(res.mean.tpr,mean(res.qt.tpr),mean(res.ca.tpr),mean(res.dum.tpr),mean(res.mix.tpr),mean(res.spls.tpr))
res.sum.td<-c(res.sum.td,sum(res.qt.td),sum(res.ca.td),sum(res.dum.td),sum(res.mix.td),sum(res.spls.td))




res.qt.res <- c(res.qt.res,max(res.qt.acc))
res.qt.res <- c(res.qt.res,mean(res.qt.acc))
res.qt.res <- c(res.qt.res,min(res.qt.acc))

res.ca.res <- c(res.ca.res,max(res.ca.acc))
res.ca.res <- c(res.ca.res,mean(res.ca.acc))
res.ca.res <- c(res.ca.res,min(res.ca.acc))

res.dum.res<- c(res.dum.res,max(res.dum.acc))
res.dum.res<- c(res.dum.res,mean(res.dum.acc))
res.dum.res<- c(res.dum.res,min(res.dum.acc))

res.mix.res <- c(res.mix.res,max(res.mix.acc))
res.mix.res <- c(res.mix.res,mean(res.mix.acc))
res.mix.res <- c(res.mix.res,min(res.mix.acc))
}
