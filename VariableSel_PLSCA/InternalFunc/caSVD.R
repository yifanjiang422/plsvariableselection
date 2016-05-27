caSVD <- function(R, masses, weights, make_data1_nominal = FALSE, make_data2_nominal = FALSE,  weights1 = NULL, weights2 = NULL,
                  DESIGN = NULL, make_design_nominal = TRUE, symmetric = TRUE, graphs = TRUE, k = 1){
#   if (nrow(DATA1) != nrow(DATA2)) {
#     stop("DATA1 and DATA2 must have the same number of rows.")
#   }
#   if (make_data1_nominal) {
#     DATA1 <- makeNominalData(DATA1)
#   }
#   if (make_data2_nominal) {
#     DATA2 <- makeNominalData(DATA2)
#   }
#   DESIGN <- texpoDesignCheck(DATA1, DESIGN, make_design_nominal)
#   DESIGN <- texpoDesignCheck(DATA2, DESIGN, make_design_nominal = FALSE)
#   DATA1 <- as.matrix(DATA1)
#   DATA2 <- as.matrix(DATA2)
  #R <- t(DATA1) %*% DATA2
  
  DATA_dimensions = dim(R)
  
  mRP <- makeRowProfiles.vs(R, weights = weights1, masses = weights2)
  pdq_results <- genPDQ(datain = mRP$deviations, M = masses, W = weights, is.mds = FALSE, decomp.approach = "svd", 
                        k = 2)
  return (list(R=R,pdq = pdq_results))
}