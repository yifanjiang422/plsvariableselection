library(MASS)
library(ggplot2)
library(reshape2)

x_axis=c(45,40,35,30,25)

#-----rho0.5-------#
qtacc05_1=c(1,1,1,1,1)
caacc05_1=c(1,1,1,1,1)
dumacc05_1=c(1,1,1,1,1)
mixacc05_1=c(0.9667,0.9375,0.9143,0.88,0.8167)

#----rho1.5-------#
qtacc15_1=c(0.953111,0.91775,0.86114,0.7773,0.6636)
caacc15_1=c(0.952889,0.9175,0.8617,0.7776,0.664)
dumacc15_1=c(0.95289,0.91775,0.8617,0.778,0.6636)
mixacc15_1=c(0.939333,0.8825,0.8248,0.7406,0.67)

#---rho2.5-------#
qtacc25_1=c(0.90489,0.81175,0.73257,0.62467,0.49)
caacc25_1=c(0.90511,0.8125,0.73314,0.6253,0.49)
dumacc25_1=c(0.90533,0.81325,0.73285,0.62633,0.4892)
mixacc25_1=c(0.90444,0.805,0.7137,0.61433,0.4952)


#-----------2ndcomponent----------------#
qtacc05_2=c(0.993,0.988,0.9833,0.9662,0.948)
caacc05_2=c(0.9943,0.9886,0.981,0.9696,0.9537)
dumacc05_2=c(0.9544,0.904,0.8517,0.7755,0.6697)
mixacc05_2=c(0.95556,0.95,0.8714,0.8083,0.74647)

#----rho1.5-------#
qtacc15_2=c(0.9407,0.90025,0.83485,0.74933,0.6298)
caacc15_2=c(0.9562,0.91987,0.8541,0.7603,0.6342)
dumacc15_2=c(0.93789,0.86675,0.78442,0.692,0.5874)
mixacc15_2=c(0.9398,0.8751,0.8075,0.721,0.6376)

#---rho2.5-------#
qtacc25_2=c(0.909778,0.822625,0.73871,0.64983,0.5386)
caacc25_2=c(0.91267,0.822,0.73214,0.6365,0.5282)
dumacc25_2=c(0.9104,0.82125,0.73185,0.63967,0.528)
mixacc25_2=c(0.9104,0.826625,0.7345,0.6528,0.5452)

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