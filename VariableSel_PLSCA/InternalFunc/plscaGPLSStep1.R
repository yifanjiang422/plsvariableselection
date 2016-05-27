plscaGPLSStep1 <- function (X,Z, ind.block.x, ind.block.z, sparsity.x, sparsity.z, 
                            epsilon, iter.max,h,masses,weights) 
{
  R<- t(X)%*%Z
  svd.R <- caSVD(R,masses,weights)
  u.tild.old <- svd.R$pdq$p[,1]
  v.tild.old <- svd.R$pdq$q[,1]
  u.tild.previous <- v.tild.previous <- 0
  iter <- 0
  while (((normv(u.tild.old - u.tild.previous) > epsilon)) & 
           (iter < iter.max)) {
    vecZV <- R %*% matrix(v.tild.old, ncol = 1)
    tab.ind <- c(0, ind.block.x, length(vecZV))
    res <- NULL
    for (i in 1:(length(ind.block.x) + 1)) {
      ji <- tab.ind[i + 1] - tab.ind[i]
      vecx <- vecZV[((tab.ind[i] + 1):tab.ind[i + 1])]
      res <- c(res, 2 * normv(vecx)/sqrt(ji))
    }
    if (sparsity.x == 0) 
      lambda.x <- 0
    else {
      lambda.x <- sort(res)[sparsity.x]
    }
    if (sparsity.z == 0) {
      lambda.z <- 0
    }
    else {
      vecZV <- t(R) %*% matrix(u.tild.old, ncol = 1)
      tab.ind <- c(0, ind.block.z, length(vecZV))
      res <- NULL
      for (i in 1:(length(ind.block.z) + 1)) {
        ji <- tab.ind[i + 1] - tab.ind[i]
        vecx <- vecZV[((tab.ind[i] + 1):tab.ind[i + 1])]
        res <- c(res, 2 * normv(vecx)/sqrt(ji))
      }
      lambda.z <- sort(res)[sparsity.z]
    }
    u.tild.new <- soft.thresholding.group(R %*% matrix(v.tild.old, 
                                                       ncol = 1), ind = ind.block.x, lambda = lambda.x)
    u.tild.new <- u.tild.new/sqrt(sum(u.tild.new^2))
    if (sparsity.z == 0) {
      v.tild.new <- t(R) %*% matrix(u.tild.old, ncol = 1)
    }
    else {
      v.tild.new <- soft.thresholding.group(t(R) %*% matrix(u.tild.old, 
                                                            ncol = 1), ind = ind.block.z, lambda = lambda.z)
    }
    v.tild.new <- v.tild.new/sqrt(sum(v.tild.new^2))
    u.tild.previous <- u.tild.old
    v.tild.previous <- v.tild.old
    u.tild.old <- u.tild.new
    v.tild.old <- v.tild.new
    #print(iter)
    iter <- iter + 1
    if(any(is.na(v.tild.new)) | any(is.na(u.tild.new))){
      iter=101
      u.tild.new<-u.tild.previous
      u.tild.new<-u.tild.previous
    }
  }
  res <- list(iter = iter, u.tild.new = u.tild.new, v.tild.new = v.tild.new)
}