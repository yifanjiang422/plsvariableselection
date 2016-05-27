library(MASS)
library(ggplot2)
library(reshape2)

x_axis=c(90,80,70,60,50)

#-----rho0.5-------#
qtacc05_1=c(1,1,1,1,1)
caacc05_1=c(1,1,1,1,1)
dumacc05_1=c(0.99989,0.99988,0.99986,0.99967,0.99931)
mixacc05_1=c(0.9928,0.99,0.9889,0.9625,0.8167)

#----rho1.5-------#
qtacc15_1=c(0.95378,0.91175,0.8487143,0.7575,0.6402)
caacc15_1=c(0.954222,0.911625,0.8401429,0.7535,0.6392)
dumacc15_1=c(0.954111,0.91075,0.8457143,0.749,0.6322)
mixacc15_1=c(0.95467,0.90425,0.841,0.75783,0.642)

#---rho2.5-------#
qtacc25_1=c(0.905444,0.81425,0.7148571,0.606333,0.4896)
caacc25_1=c(0.905444,0.8065,0.7108571,0.6036667,0.4802)
dumacc25_1=c(0.9052222,0.812625,0.710743,0.6015,0.4688)
mixacc25_1=c(0.90633,0.81125,0.70986,0.6225,0.4836)


#-----------2ndcomponent----------------#
qtacc05_2=c(0.9949444,0.989875,0.9826429,0.972,0.96967)
caacc05_2=c(0.9938333,0.9849375,0.9639286,0.95333,0.9369)
dumacc05_2=c(0.9627222,0.9195625,0.8527857,0.7945,0.6735)
mixacc05_2=c(0.9695,0.9268,0.8681,0.8126,0.7094)

#----rho1.5-------#
qtacc15_2=c(0.934222,0.8961875,0.8374286,0.7315,0.6066)
caacc15_2=c(0.947778,0.88175,0.8130714,0.7256667,0.6025)
dumacc15_2=c(0.931444,0.8654375,0.7851429,0.70075,0.5866)
mixacc15_2=c(0.9467,0.8956,0.8144,0.7199,0.5937)

#---rho2.5-------#
qtacc25_2=c(0.9075,0.8196875,0.7262143,0.6254167,0.5187)
caacc25_2=c(0.9062778,0.8189375,0.7220714,0.63025,0.5135)
dumacc25_2=c(0.9000162,0.818375,0.7211429,0.6211667,0.5104)
mixacc25_2=c(0.9112,0.8271,0.73114,0.63325,0.5146)

#------plotrho05comp1------#
plotdata<-data.frame(
  splsos = qtacc05_1,
  gplsca = caacc05_1,
  gplsdum = dumacc05_1,
  mix=mixacc05_1,
  xvalue = x_axis)
plotdatalong <- melt(plotdata,id="xvalue")
rho05results <- ggplot(plotdatalong,
                     aes(x=xvalue,y=value,colour=variable,shape=variable))+geom_point()+geom_line()+scale_shape_manual(values=c(1,2,3,4))+labs(title="Variable Selection Simulation Performance (rho=0.5)",x="Sparsity",y="Accuracy") 
rho05results=rho05results+ylim(0.4,1)
print(rho05results)

#------plotrho15comp1------#
plotdata<-data.frame(
  splsos = qtacc15_1,
  gplsca = caacc15_1,
  gplsdum = dumacc15_1,
  mix=mixacc15_1,
  xvalue = x_axis)
plotdatalong <- melt(plotdata,id="xvalue")
rho15results <- ggplot(plotdatalong,
                       aes(x=xvalue,y=value,colour=variable,shape=variable))+geom_point()+geom_line()+scale_shape_manual(values=c(1,2,3,4))+labs(title="Variable Selection Simulation Performance (rho=1.5)",x="Sparsity",y="Accuracy") 
rho15results=rho15results+ylim(0.4,1)
print(rho15results)

#------plotrho25comp1------#
plotdata<-data.frame(
  splsos = qtacc25_1,
  gplsca = caacc25_1,
  gplsdum = dumacc25_1,
  mix=mixacc25_1,
  xvalue = x_axis)
plotdatalong <- melt(plotdata,id="xvalue")
rho25results <- ggplot(plotdatalong,
                       aes(x=xvalue,y=value,colour=variable,shape=variable))+geom_point()+geom_line()+scale_shape_manual(values=c(1,2,3,4))+labs(title="Variable Selection Simulation Performance (rho=2.5)",x="Sparsity",y="Accuracy") 
rho25results=rho25results+ylim(0.4,1)
print(rho25results)

#------plotrho05comp2------#
plotdata<-data.frame(
  splsos = qtacc05_2,
  gplsca = caacc05_2,
  gplsdum = dumacc05_2,
  mix=mixacc05_2,
  xvalue = x_axis)
plotdatalong <- melt(plotdata,id="xvalue")
rho05results <- ggplot(plotdatalong,
                       aes(x=xvalue,y=value,colour=variable,shape=variable))+geom_point()+geom_line()+scale_shape_manual(values=c(1,2,3,4))+labs(title="Variable Selection Simulation Performance (rho=0.5, ncomp=2)",x="Sparsity",y="Accuracy") 
rho05results=rho05results+ylim(0.4,1)
print(rho05results)

#------plotrho15comp2------#
plotdata<-data.frame(
  splsos = qtacc15_2,
  gplsca = caacc15_2,
  gplsdum = dumacc15_2,
  mix=mixacc15_2,
  xvalue = x_axis)
plotdatalong <- melt(plotdata,id="xvalue")
rho15results <- ggplot(plotdatalong,
                       aes(x=xvalue,y=value,colour=variable,shape=variable))+geom_point()+geom_line()+scale_shape_manual(values=c(1,2,3,4))+labs(title="Variable Selection Simulation Performance (rho=1.5,ncomp=2)",x="Sparsity",y="Accuracy") 
rho15results=rho15results+ylim(0.4,1)
print(rho15results)

#------plotrho25comp2------#
plotdata<-data.frame(
  splsos = qtacc25_2,
  gplsca = caacc25_2,
  gplsdum = dumacc25_2,
  mix=mixacc25_2,
  xvalue = x_axis)
plotdatalong <- melt(plotdata,id="xvalue")
rho25results <- ggplot(plotdatalong,
                       aes(x=xvalue,y=value,colour=variable,shape=variable))+geom_point()+geom_line()+scale_shape_manual(values=c(1,2,3,4))+labs(title="Variable Selection Simulation Performance (rho=2.5,ncomp=2)",x="Sparsity",y="Accuracy") 
rho25results=rho25results+ylim(0.4,1)
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
