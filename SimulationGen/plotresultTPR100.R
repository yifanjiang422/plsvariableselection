library(MASS)
library(ggplot2)
library(reshape2)

x_axis=c(90,70,50,30,10)

#-----rho0.5-------#
qtacc05_1=c(1,1,1,1,1)
caacc05_1=c(1,1,1,1,1)
dumacc05_1=c(1,1,1,1,1)
mixacc05_1=c(0.88,0.911,0.9292,0.9557,0.9827)
splsacc05_1=c(1,1,1,1,1)

#----rho1.5-------#
qtacc15_1=c(0.596,0.6463,0.6322,0.6841,0.8905)
caacc15_1=c(0.597,0.6473,0.6304,0.6831,0.89044)
dumacc15_1=c(0.597,0.6746,0.6294,0.6835,0.89167)
mixacc15_1=c(0.477,0.575,0.6432,0.7218,0.897)
splsacc15_1=c(0.796,0.72767,0.6812,0.739,0.9034)
#---rho2.5-------#
qtacc25_1=c(0.168,0.3467,0.4796,0.682,0.7464)
caacc25_1=c(0.168,0.3476,0.478,0.682,0.7429)
dumacc25_1=c(0.168,0.34833,0.477,0.6789,0.7395)
mixacc25_1=c(0.116,0.3247,0.4866,0.6913,0.8717)
splsacc25_1=c(0.222,0.4183,0.5026,0.7028,0.8993)

#-----------2ndcomponent----------------#
qtacc05_2=c(0.9625,0.9488,0.6471,0.7375,0.9054)
caacc05_2=c(0.971,0.9665,0.6446,0.729413,0.9034)
dumacc05_2=c(0.609,0.6492,0.5491,0.7003,0.9006)
mixacc05_2=c(0.7435,0.7543,0.6138,0.746,0.9033)
splsacc05_2=c(0.9895,0.988,0.6728,0.761,0.9069)
#----rho1.5-------#
qtacc15_2=c(0.5275,0.5818,0.512,0.69557,0.8947)
caacc15_2=c(0.649,0.6243,0.5172,0.68864,0.8972)
dumacc15_2=c(0.41,0.4941,0.516,0.693,0.898)
mixacc15_2=c(0.4155,0.527,0.5138,0.7115,0.8968)
splsacc15_2=c(0.7855,0.6985,0.5265,0.73521,0.8994)

#---rho2.5-------#
qtacc25_2=c(0.1905,0.3425,0.5103,0.708,0.8367)
caacc25_2=c(0.192,0.3557,0.5148,0.7135,0.8451)
dumacc25_2=c(0.181,0.3517,0.5153,0.7091,0.857)
mixacc25_2=c(0.164,0.341,0.5103,0.694,0.8914)
splsacc25_2=c(0.228,0.4215,0.5206,0.7246,0.8978)

#------plotrho05comp1------#
plotdata<-data.frame(
  splsos = qtacc05_1,
  gplsca = caacc05_1,
  gplsdum = dumacc05_1,
  mix=mixacc05_1,
  spls=splsacc05_1,
  xvalue = x_axis)
plotdatalong <- melt(plotdata,id="xvalue")
rho05results <- ggplot(plotdatalong,
                     aes(x=xvalue,y=value,colour=variable,shape=variable))+geom_point()+geom_line()+scale_shape_manual(values=c(1,2,3,4,5))+labs(title="Variable Selection Simulation Performance (sigma=0.5, ncomp=1)",x="Sparsity",y="Accuracy") 
rho05results=rho05results+ylim(0,1)
print(rho05results)

#------plotrho15comp1------#
plotdata<-data.frame(
  splsos = qtacc15_1,
  gplsca = caacc15_1,
  gplsdum = dumacc15_1,
  mix=mixacc15_1,
  spls=splsacc15_1,
  xvalue = x_axis)
plotdatalong <- melt(plotdata,id="xvalue")
rho15results <- ggplot(plotdatalong,
                       aes(x=xvalue,y=value,colour=variable,shape=variable))+geom_point()+geom_line()+scale_shape_manual(values=c(1,2,3,4,5))+labs(title="Variable Selection Simulation Performance (rho=1.5, ncomp=1)",x="Sparsity",y="Accuracy") 
rho15results=rho15results+ylim(0,1)
print(rho15results)

#------plotrho25comp1------#
plotdata<-data.frame(
  splsos = qtacc25_1,
  gplsca = caacc25_1,
  gplsdum = dumacc25_1,
  mix=mixacc25_1,
  spls=splsacc25_1,
  xvalue = x_axis)
plotdatalong <- melt(plotdata,id="xvalue")
rho25results <- ggplot(plotdatalong,
                       aes(x=xvalue,y=value,colour=variable,shape=variable))+geom_point()+geom_line()+scale_shape_manual(values=c(1,2,3,4,5))+labs(title="Variable Selection Simulation Performance (rho=2.5,ncomp=1)",x="Sparsity",y="Accuracy") 
rho25results=rho25results+ylim(0,1)
print(rho25results)

#------plotrho05comp2------#
plotdata<-data.frame(
  splsos = qtacc05_2,
  gplsca = caacc05_2,
  gplsdum = dumacc05_2,
  mix=mixacc05_2,
  spls=splsacc05_2,
  xvalue = x_axis)
plotdatalong <- melt(plotdata,id="xvalue")
rho05results <- ggplot(plotdatalong,
                       aes(x=xvalue,y=value,colour=variable,shape=variable))+geom_point()+geom_line()+scale_shape_manual(values=c(1,2,3,4,5))+labs(title="Variable Selection Simulation Performance (rho=0.5, ncomp=2)",x="Sparsity",y="Accuracy") 
rho05results=rho05results+ylim(0,1)
print(rho05results)

#------plotrho15comp2------#
plotdata<-data.frame(
  splsos = qtacc15_2,
  gplsca = caacc15_2,
  gplsdum = dumacc15_2,
  mix=mixacc15_2,
  spls=splsacc15_1,
  xvalue = x_axis)
plotdatalong <- melt(plotdata,id="xvalue")
rho15results <- ggplot(plotdatalong,
                       aes(x=xvalue,y=value,colour=variable,shape=variable))+geom_point()+geom_line()+scale_shape_manual(values=c(1,2,3,4,5))+labs(title="Variable Selection Simulation Performance (rho=1.5,ncomp=2)",x="Sparsity",y="Accuracy") 
rho15results=rho15results+ylim(0,1)
print(rho15results)

#------plotrho25comp2------#
plotdata<-data.frame(
  splsos = qtacc25_2,
  gplsca = caacc25_2,
  gplsdum = dumacc25_2,
  mix=mixacc25_2,
  spls=splsacc25_1,
  xvalue = x_axis)
plotdatalong <- melt(plotdata,id="xvalue")
rho25results <- ggplot(plotdatalong,
                       aes(x=xvalue,y=value,colour=variable,shape=variable))+geom_point()+geom_line()+scale_shape_manual(values=c(1,2,3,4,5))+labs(title="Variable Selection Simulation Performance (rho=2.5,ncomp=2)",x="Sparsity",y="Accuracy") 
rho25results=rho25results+ylim(0,1)
print(rho25results)
# #------plotrho05comp2------#
# plotdata<-data.frame(
#   splsos_rho05_1 = qtacc05_1,
#   splsos_rho15_1 = qtacc15_1,
#   splsos_rho25_1 = qtacc25_1,
#   splsos_rho05_2 = qtacc05_2,
#   splsos_rho15_2 = qtacc15_2,
#   splsos_rho25_2 = qtacc25_2,
#   xvalue = x_axis)
# plotdatalong <- melt(plotdata,id="xvalue")
# rhoqtresults <- ggplot(plotdatalong,
#                        aes(x=xvalue,y=value,colour=variable,shape=variable))+geom_point()+geom_line()+scale_shape_manual(values=c(1,2,3,4,5,6))+labs(title="SPLSOS Performance witth Two LVs",x="Sparsity",y="Accuracy") 
# print(rhoqtresults)
# 
# #------plotrho15comp2------#
# plotdata<-data.frame(
#   gplsca_rho05_1 = caacc05_1,
#   gplsca_rho15_1 = caacc15_1,
#   gplsca_rho25_1 = caacc25_1,
#   gplsca_rho05_2 = caacc05_2,
#   gplsca_rho15_2 = caacc15_2,
#   gplsca_rho25_2 = caacc25_2,
#   xvalue = x_axis)
# plotdatalong <- melt(plotdata,id="xvalue")
# rhocaresults <- ggplot(plotdatalong,
#                        aes(x=xvalue,y=value,colour=variable,shape=variable))+geom_point()+geom_line()+scale_shape_manual(values=c(1,2,3,4,5,6))+labs(title="GPLSCA Performance with Two LVs",x="Sparsity",y="Accuracy") 
# print(rhocaresults)
# 
# #------plotrho25comp2------#
# plotdata<-data.frame(
#   gplsdum_rho05_1 = dumacc05_1,
#   gplsdum_rho15_1 = dumacc15_1,
#   gplsdum_rho25_1 = dumacc25_1,
#   gplsdum_rho05_2 = dumacc05_2,
#   gplsdum_rho15_2 = dumacc15_2,
#   gplsdum_rho25_2 = dumacc25_2,
#   xvalue = x_axis)
# plotdatalong <- melt(plotdata,id="xvalue")
# rhodumresults <- ggplot(plotdatalong,
#                        aes(x=xvalue,y=value,colour=variable,shape=variable))+geom_point()+geom_line()+scale_shape_manual(values=c(1,2,3,4,5,6))+labs(title="GPLSDUM Selection Simulation Performance with Two LVs",x="Sparsity",y="Accuracy") 
# print(rhodumresults)
