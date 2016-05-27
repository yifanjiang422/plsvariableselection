plscaGPLS <- function (X, Z, ncomp, mode = "regression", max.iter = 100, tol = 1e-06, 
          keepX, keepZ = NULL, ind.block.x, ind.block.z = NULL, masses, weights) 
{
  X <- as.matrix(X)
  Z <- as.matrix(Z)
 
  q <- ncol(Z)
  p <- ncol(X)
  n <- nrow(X)
  #Scale X and Z
  X.s <- (X - matrix(rep(colSums(X)/n,n),n,p,byrow=TRUE))/p
  Z.s <- (Z - matrix(rep(colSums(Z)/n,n),n,q,byrow=TRUE))/q
  #End of scale
  X.names = dimnames(X)[[2]]
  if (is.null(X.names)) 
    X.names = paste("X", 1:p, sep = "")
  if (dim(Z)[2] == 1) 
    Z.names = "Z"
  else {
    Z.names = dimnames(Z)[[2]]
    if (is.null(Z.names)) 
      Z.names = paste("Z", 1:q, sep = "")
  }
  ind.names = dimnames(X)[[1]]
  if (is.null(ind.names)) {
    ind.names = dimnames(Z)[[1]]
    rownames(X) = ind.names
  }
  if (is.null(ind.names)) {
    ind.names = 1:n
    rownames(X) = rownames(Z) = ind.names
  }
  mat.c <- matrix(nrow = p, ncol = ncomp)
  mat.d <- matrix(nrow = q, ncol = ncomp)
  mat.e <- matrix(nrow = q, ncol = ncomp)
  #X.s <- scale(X)
  #Z.s <- scale(Z)
  sparsity.x <- length(ind.block.x) + 1 - keepX
  if (is.null(ind.block.z)) {
    sparsity.z <- rep(0, ncomp)
  }
  else {
    if (is.null(keepZ)) 
      keepZ <- rep(length(ind.block.z) + 1, ncomp)
    sparsity.z <- length(ind.block.z) + 1 - keepZ
  }
  #R <- t(X) %*% Z
  res.load <- plscaGPLSStep1(X,Z, ind.block.x = ind.block.x, 
                                        ind.block.z = ind.block.z, sparsity.x = sparsity.x[1], 
                                        sparsity.z = sparsity.z[1], epsilon = tol, iter.max = max.iter,h=1,
                             masses,weights)
  res.deflat <- plscaGPLSStep2(X=X.s,Z=Z.s, res.load$u.tild.new, 
                           res.load$v.tild.new, mode = mode)
  mat.c[, 1] <- res.deflat$c
  iter <- res.load$iter
  if (mode == "regression") 
    mat.d[, 1] <- res.deflat$d
  else mat.e[, 1] <- res.deflat$e
  load.u <- res.load$u.tild.new
  load.v <- res.load$v.tild.new
  if (ncomp > 1) {
    for (h in 2:ncomp) {
      res.load <- plscaGPLSStep1(X = res.deflat$X.h, 
                                            Z = res.deflat$Z.h, ind.block.x = ind.block.x, 
                                            ind.block.z = ind.block.z, sparsity.x = sparsity.x[h], 
                                            sparsity.z = sparsity.z[h], epsilon = tol, iter.max = max.iter,h=h,masses,weights)
      load.u <- cbind(load.u, res.load$u.tild.new)
      load.v <- cbind(load.v, res.load$v.tild.new)
      res.deflat <- plscaGPLSStep2(X = res.deflat$X.h, Z = res.deflat$Z.h, 
                               res.load$u.tild.new, res.load$v.tild.new, mode = mode)
      mat.c[, h] <- res.deflat$c
      if (mode == "regression") 
        mat.d[, h] <- res.deflat$d
      else mat.e[, h] <- res.deflat$e
      iter <- c(iter, res.load$iter)
    }
  }
  else {
    load.u <- matrix(load.u, ncol = 1)
    load.v <- matrix(load.v, ncol = 1)
  }
  mat.t <- X.s %*% load.u
  mat.u <- Z.s %*% load.v
  cl = match.call()
  if (is.null(keepZ)) {
    if (is.null(ind.block.z)) 
      keepZ <- rep(ncol(Z), ncomp)
    else keepZ <- rep(length(ind.block.z) + 1, ncomp)
  }
  result <- list(call = cl, X = X, Z = Z, ncomp = ncomp, 
                 mode = mode, keepX = keepX, keepZ = keepZ, mat.c = mat.c, 
                 mat.d = mat.d, mat.e = mat.e, loadings = list(X = load.u, 
                                                               Z = load.v), variates = list(X = mat.t, Z = mat.u), 
                 names = list(X = X.names, Z = Z.names, indiv = ind.names), 
                 tol = tol, max.iter = max.iter, iter = iter, ind.block.x = ind.block.x, 
                 ind.block.z = ind.block.z)
  class(result) = c("gPLS", "sPLS", "spls", "pls")
  return(invisible(result))
}