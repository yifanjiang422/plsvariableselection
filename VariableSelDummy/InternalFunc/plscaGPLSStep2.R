plscaGPLSStep2 <- function (X,Z, u.tild.new, v.tild.new, mode) 
{
  
  xi.h <- X %*% matrix(u.tild.new, ncol = 1)/((normv(u.tild.new))^2)
  w.h <- Z %*% matrix(v.tild.new, ncol = 1)/((normv(v.tild.new))^2)
  c.h <- t(X) %*% matrix(xi.h, ncol = 1)/((normv(xi.h))^2)
  d.rh <- t(Z) %*% matrix(xi.h, ncol = 1)/(sum(xi.h * xi.h))
  d.h <- t(Z) %*% matrix(w.h, ncol = 1)/(sum(w.h * w.h))
  X.h <- X - xi.h %*% t(c.h)
  if (mode == "regression") 
    Z.h <- Z - xi.h %*% t(d.rh)
  else Z.h <- Z - w.h %*% t(d.h)
  res <- list(X.h = X.h, Z.h = Z.h, c = c.h, d = d.rh, e = d.h)
  return(res)
}