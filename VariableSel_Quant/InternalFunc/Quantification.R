osquant <- function(RawX,RawZ,ncomp,quantscale,data1.qnt=NULL,data2.qnt=NULL){
  #RawData.Matrix<-as.matrix(Data.Raw[,1:No.Var])           #Load matrix X as the variable matrix
  #RawX <- beer.tasting.notes$data[,1:8]
  #RawZ <- as.matrix(beer.tasting.notes$data[,11:13])
  #RawX <- data.matrix(beer.tasting.notes$city.design)       #Define X matrix                
  #RawZ <- data.matrix(beer.tasting.notes$pale.sour.style)   #Define Z matrix
  
  if(is.null(data1.qnt)){
    p<-ncol(RawX)
    }#Define number of cat var in X
  else{
    p <- ncol(RawX)+ncol(data1.qnt)
    data1.qnt<-scale(data1.qnt)*quantscale
  }
  
  if (is.null(data2.qnt)){
    q=ncol(RawZ)
  }else{
    q <- ncol(RawZ)+ncol(data2.qnt)                                           #Define number of cat var in Z
    data2.qnt<-scale(data2.qnt)*quantscale
  }

  colnameRawX<-colnames(RawX)                               #store column variable name X
  colnameRawZ<-colnames(RawZ)                               #store column variable name Z
  
 
  #Non-metric Variable quantification#############
  #ncomp<-3                                                  #No of components = no of lcs to be confirmed with GPLS algorithm = 5 (assumes)
  noiteration<-0                                            #initialized constant                        
  n<-nrow(RawX)                                             #Number of observations
  
  #Omega<-matrix(c(rep(0,n*ncomp)),n,ncomp)                                   #Y Score, component in the Y space 
  Omega<-matrix(c(rep(c(1,rep(0,(n-1))),ncomp)),n,ncomp)
#   if(is.null(data2.qnt)){
#   }else{
#     Omega[,1]=data2.qnt[,1]
#   }
  U<-matrix(c(rep(0,p*ncomp)),p,ncomp)                                       #X weight X loading
  Xi<-matrix(c(rep(0,n*ncomp)),n,ncomp)                                      #X score, Latend Component as a linear component of X variables
  betaXi<-matrix(c(rep(0,n*ncomp)),n,ncomp)
  V<-matrix(c(rep(0,q*ncomp)),q,ncomp)                                       #Y Weight                            
  correc<- quantscale
  #beta<-matrix(c(rep(0,ncomp)),ncomp,1)
  #Pmat<-matrix(c(rep(0,p*ncomp)),p,ncomp)
  #Xarray<-array(,c(n,p,ncomp))
  #Zarray<-array(,c(n,q,ncomp))
  
  i=0
  repeat{
    OmegaStart<-Omega[,i+1]
    if(i==0){
      RawX.quant <- GPLSos(Omega[,i+1],RawX)
      colnames(RawX.quant)<-colnameRawX
      #QuantX <- scale(RawX.quant)*correc
      QuantX_org <-scale(RawX.quant)*correc
    }
    if(is.null(data1.qnt)){
      QuantX <- scale(RawX.quant)*correc
    }#Define number of cat var in X
    else{
      QuantX <- cbind(scale(RawX.quant)*correc,data1.qnt)
    }
    U[,i+1]<-as.matrix((t(QuantX)%*%Omega[,i+1])/as.numeric(t(Omega[,i+1])%*%Omega[,i+1]))
    U[,i+1]<-U[,i+1]/sqrt(as.numeric(t(U[,i+1])%*%U[,i+1]))
    Xi[,i+1]<-(QuantX%*%U[,i+1])/as.numeric(t(U[,i+1])%*%U[,i+1])
    if(i==0){
      RawZ.quant <- GPLSos(Xi[,i+1],RawZ)
      colnames(RawZ.quant)<-colnameRawZ
      #QuantZ<-scale(RawZ.quant)*correc
      QuantZ_org<-scale(RawZ.quant)*correc
    }
    if (is.null(data2.qnt)){
      QuantZ <- scale(RawZ.quant)*correc
    }else{
      QuantZ <- cbind(data2.qnt,scale(RawZ.quant)*correc)
    }
    QuantZ <- cbind(scale(RawZ.quant)*correc,data2.qnt)
    V[,i+1]<-(t(QuantZ)%*%Xi[,i+1])/as.numeric(norm(t(Xi[,i+1])%*%Xi[,i+1]))
    V[,i+1]<-V[,i+1]/sqrt(as.numeric(t(V[,i+1])%*%V[,i+1]))
    Omega[,i+1]<-(QuantZ%*%V[,i+1])/as.numeric(t(V[,i+1])%*%V[,i+1])
    conv<-max(abs(OmegaStart-Omega[,i+1]))
    #print("Conv:");print(conv)
    noiteration<-noiteration+1
    if(conv<0.0000001 | noiteration>149){
      break
    }
  }
  
  return(res.quant <- list(QuantX=QuantX, QuantZ=QuantZ, Omega=Omega, Xi=Xi, U=U, V=V,noiteration=noiteration))
}




