makeRowProfiles.vs <- function (X, weights = NULL, masses = NULL, hellinger = FALSE,h) 
{
  X_dimensions <- dim(X)
  colTotal <- colSums(X)
  rowTotal <- rowSums(X)
  rowTotal[which(rowTotal==0)]<-1
  grandTotal <- sum(X)
  if (hellinger) {
    return(hellingerNorm(X, X_dimensions, colTotal, rowTotal, 
                         grandTotal, weights, masses))
  }
  else {
    return(caNormvs(X, X_dimensions, colTotal, rowTotal, grandTotal, 
                  weights, masses,h))
  }
}