Quant_step1.spls.sparsity<- function (X, Y, sparsity.x, sparsity.y, epsilon, iter.max, Omega,Xi,U,V) {
  n <- dim(X)[1]
  Z <- t(X) %*% Y
  #svd.Z <- svd(Z, nu = 1, nv = 1)
  u.tild.old <- U[,1]
  v.tild.old <- V[,1]
  u.tild.previous <- v.tild.previous <- 0
  iter <- 0
  while (((normv(u.tild.old - u.tild.previous)/normv(u.tild.old)) > 
          epsilon) & (iter < iter.max)) {
    if (sparsity.x == 0) {
      lambda.x <- 0
    }
    else {
      lambda.x <- sort(abs(Z %*% matrix(v.tild.old, ncol = 1)))[sparsity.x]
    }
    u.tild.new <- soft.thresholding(Z %*% matrix(v.tild.old, 
                                                 ncol = 1), lambda = lambda.x)
    u.tild.new <- u.tild.new/sqrt(sum(u.tild.new^2))
    if (sparsity.y == 0) 
      lambda.y <- 0
    else lambda.y <- sort(abs(t(Z) %*% matrix(u.tild.old, 
                                              ncol = 1)))[sparsity.y]
    v.tild.new <- soft.thresholding(t(Z) %*% matrix(u.tild.old, 
                                                    ncol = 1), lambda = lambda.y)
    v.tild.new <- v.tild.new/sqrt(sum(v.tild.new^2))
    u.tild.previous <- u.tild.old
    v.tild.previous <- v.tild.old
    u.tild.old <- u.tild.new
    v.tild.old <- v.tild.new
    iter <- iter + 1
  }
  res <- list(iter = iter, u.tild.new = u.tild.new, v.tild.new = v.tild.new)
}
  