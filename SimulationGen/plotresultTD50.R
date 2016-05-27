library(MASS)
library(ggplot2)
library(reshape2)

x_axis=c(45,30,20,10,5)

#-----rho0.5-------#
qtacc05_1=c(0,0,0,0,0)
caacc05_1=c(0,0,0,0,0)
dumacc05_1=c(0,0,0,0,0)
mixacc05_1=c(78,208,184,124,72)
splsacc05_1=c(0,0,0,0,0)

#----rho1.5-------#
qtacc15_1=c(388,1333,1824,1622,1094)
caacc15_1=c(1164,3957,5457,4872,3273)
dumacc15_1=c(1170,3966,5472,4896,3273)
mixacc15_1=c(496,1466,1773,1454,885)
splsacc15_1=c(224,1072,1678,1466,844)
#---rho2.5-------#
qtacc25_1=c(846,2234,2300,1876,1084)
caacc25_1=c(2532,6696,6897,5598,3255)
dumacc25_1=c(2532,6702,6900,5640,3279)
mixacc25_1=c(872,2238,2222,1630,945)
splsacc25_1=c(782,2008,2182,1505,894)

#-----------2ndcomponent----------------#
qtacc05_2=c(60,266,3804,2916,1688)
caacc05_2=c(150,816,11748,9126,5046)
dumacc05_2=c(2400,8184,13446,9396,5304)
mixacc05_2=c(502,2168,3768,2798,1740)
splsacc05_2=c(54,118,3478,2746,1656)
#----rho1.5-------#
qtacc15_2=c(1012,3055,4461,3071,1918)
caacc15_2=c(2346,8511,13401,9012,5655)
dumacc15_2=c(3414,10889,13575,9099,5628)
mixacc15_2=c(1170,3196,4573,2915,1827)
splsacc15_2=c(386,2386,4526,2920,1774)

#---rho2.5-------#
qtacc25_2=c(1592,4310,4423,3172,1972)
caacc25_2=c(4752,12789,13296,9879,6057)
dumacc25_2=c(4872,12846,13503,9609,5961)
mixacc25_2=c(1590,4296,4586,2975,1872)
splsacc25_2=c(1508,4092,4420,2990,1842)

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
rho05results=rho05results+ylim(0,14000)
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
rho15results=rho15results+ylim(0,14000)
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
rho25results=rho25results+ylim(0,14000)
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
rho05results=rho05results+ylim(0,14000)
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
rho15results=rho15results+ylim(0,14000)
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
rho25results=rho25results+ylim(0,14000)
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
