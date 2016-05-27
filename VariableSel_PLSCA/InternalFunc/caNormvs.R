caNormvs <- function (X, X_dimensions, colTotal, rowTotal, grandTotal, weights = NULL, 
          masses = NULL,h) 
{
  rowCenter = colTotal/grandTotal
  if (is.null(masses)) {
    masses = rowTotal/grandTotal
  }
  if (is.null(weights)) {
    weights = rowCenter^-1
  }
 # X<-ifelse(X==0,0.001,X)
  rowProfiles <- rowNormsvs(X, type = "ca",rowTotal)
  
  deviations <- rowProfiles - matrix(rowCenter, X_dimensions[1], 
                                     X_dimensions[2], byrow = TRUE)
  return(list(rowCenter = rowCenter, masses = masses, weights = weights, 
              rowProfiles = rowProfiles, deviations = deviations))
}