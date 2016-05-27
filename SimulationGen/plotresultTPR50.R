library(MASS)
library(ggplot2)
library(reshape2)

x_axis=c(45,30,20,10,5)

#-----rho0.5-------#
qtacc05_1=c(1,1,1,1,1)
caacc05_1=c(1,1,1,1,1)
dumacc05_1=c(1,1,1,1,1)
mixacc05_1=c(0.922,0.948,0.9693,0.9845,0.992)
splsacc05_1=c(1,1,1,1,1)

#----rho1.5-------#
qtacc15_1=c(0.612,0.665,0.688,0.76025,0.8473)
caacc15_1=c(0.612,0.668,0.6893,0.75875,0.848)
dumacc15_1=c(0.61,0.6675,0.688333,0.7562,0.848)
mixacc15_1=c(0.504,0.6325,0.69867,0.80875,0.89733)
splsacc15_1=c(0.776,0.732,0.7203,0.81675,0.90622)
#---rho2.5-------#
qtacc25_1=c(0.154,0.427,0.5743,0.69925,0.8556)
caacc25_1=c(0.156,0.4265,0.5723,0.7025,0.85533)
dumacc25_1=c(0.156,0.4225,0.5737,0.697,0.85267)
mixacc25_1=c(0.128,0.43,0.61,0.78075,0.89044)
splsacc25_1=c(0.218,0.498,0.633,0.8115,0.90067)

#-----------2ndcomponent----------------#
qtacc05_2=c(0.97,0.9688,0.683,0.8178,0.9062)
caacc05_2=c(0.975,0.966,0.67367,0.80987,0.9065)
dumacc05_2=c(0.6,0.659,0.6265,0.8043,0.9017)
mixacc05_2=c(0.749,0.729,0.686,0.825,0.903)
splsacc05_2=c(0.973,0.9852,0.7102,0.8284,0.908)
#----rho1.5-------#
qtacc15_2=c(0.494,0.61725,0.621,0.79125,0.884667)
caacc15_2=c(0.609,0.6445,0.62067,0.786875,0.88778)
dumacc15_2=c(0.431,0.5455,0.617,0.789,0.89167)
mixacc15_2=c(0.415,0.59975,0.6155,0.80525,0.8952)
splsacc15_2=c(0.807,0.70175,0.6228,0.8175,0.90144)

#---rho2.5-------#
qtacc25_2=c(0.204,0.44875,0.60767,0.739625,0.87933)
caacc25_2=c(0.208,0.4565,0.60967,0.731875,0.876)
dumacc25_2=c(0.188,0.452,0.6025,0.743625,0.88011)
mixacc25_2=c(0.205,0.44875,0.59817,0.78475,0.895778)
splsacc25_2=c(0.246,0.4885,0.63167,0.813125,0.89767)

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
