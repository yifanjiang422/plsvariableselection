quant.step2.spls <- function (X, Y, u.tild.new, v.tild.new, Xi, Omega, mode) 
{
  xi.h <- Xi
  w.h <- Omega
  c.h <- t(X) %*% matrix(xi.h, ncol = 1)/((normv(xi.h))^2)
  d.rh <- t(Y) %*% matrix(xi.h, ncol = 1)/(sum(xi.h * xi.h))
  d.h <- t(Y) %*% matrix(w.h, ncol = 1)/(sum(w.h * w.h))
  X.h <- X - xi.h %*% t(c.h)
  if (mode == "regression") 
    Y.h <- Y - xi.h %*% t(d.rh)
  else Y.h <- Y - w.h %*% t(d.h)
  res <- list(X.h = X.h, Y.h = Y.h, c = c.h, d = d.rh, e = d.h)
  return(res)
}