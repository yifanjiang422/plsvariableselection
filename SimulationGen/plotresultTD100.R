library(MASS)
library(ggplot2)
library(reshape2)

x_axis=c(90,70,50,30,10)

#-----rho0.5-------#
qtacc05_1=c(0,0,0,0,0)
caacc05_1=c(0,0,0,0,0)
dumacc05_1=c(0,0,0,0,0)
mixacc05_1=c(240,534,708,310,72)
splsacc05_1=c(0,0,0,0,0)

#----rho1.5-------#
qtacc15_1=c(808,2099,3404,3861,1884)
caacc15_1=c(806,2096,3409,3867,1884)
dumacc15_1=c(806,2095,3409,3866,1876)
mixacc15_1=c(1046,2519,3235,3606,1832)
splsacc15_1=c(408,1634,3188,3654,1738)
#---rho2.5-------#
qtacc25_1=c(1664,3782,5002,4309,2995)
caacc25_1=c(1664,3778,5003,4306,3023)
dumacc25_1=c(1664,3774,5000,4314,3049)
mixacc25_1=c(1768,3929,4990,4258,2019)
splsacc25_1=c(1556,3490,4974,4160,1812)

#-----------2ndcomponent----------------#
qtacc05_2=c(150,614,7058,7350,3406)
caacc05_2=c(116,402,7108,7584,3478)
dumacc05_2=c(1564,4210,9018,8392,3578)
mixacc05_2=c(1026,2948,7724,7112,3480)
splsacc05_2=c(42,144,6544,6688,3350)
#----rho1.5-------#
qtacc15_2=c(1890,4977,9513,7559,3614)
caacc15_2=c(1404,4474,9402,7815,3565)
dumacc15_2=c(2360,6032,9425,7714,3552)
mixacc15_2=c(2338,5627,9411,7431,3649)
splsacc15_2=c(858,3618,9470,7414,3622)

#---rho2.5-------#
qtacc25_2=c(3238,7634,9431,7913,4484)
caacc25_2=c(3230,7526,9420,7784,4413)
dumacc25_2=c(3275,7489,9387,7875,4260)
mixacc25_2=c(3344,7704,9480,8468,3796)
splsacc25_2=c(3088,6942,9588,7710,3680)

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
rho05results=rho05results+ylim(0,10000)
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
rho15results=rho15results+ylim(0,10000)
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
rho25results=rho25results+ylim(0,10000)
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
rho05results=rho05results+ylim(0,10000)
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
rho15results=rho15results+ylim(0,10000)
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
rho25results=rho25results+ylim(0,10000)
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
